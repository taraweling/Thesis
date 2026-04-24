import csv # do i need this?
from pyvis.network import Network as net
import networkx as nx
from networkx.algorithms import bipartite as nx_bipartite
import matplotlib.pyplot as plt
import numpy as np
import os
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import math

try:
    from scipy.stats import mannwhitneyu
except Exception:
    mannwhitneyu = None

# Helpers

def _extract_deg_genes(degset):
    """
    Normalize DEG inputs into a set of gene IDs.

    Accepts:
      - set/iterable of gene strings
      - list of DEG records (gene in index 0)
      - dict keyed by gene
    """
    if degset is None:
        return set()

    # dict: keys are genes
    if isinstance(degset, dict):
        return set(degset.keys())

    # set/tuple/list of strings OR records
    try:
        degset_iter = list(degset)
    except TypeError:
        return set()

    if not degset_iter:
        return set()

    first = degset_iter[0]
    if isinstance(first, str):
        return set(degset_iter)

    if isinstance(first, (list, tuple)) and first:
        return {rec[0] for rec in degset_iter if rec}

    return set()

def _edge_log2fc_vals(edge, gene, gene_to_disorders):
    """
    Return log2fc values for a given edge. Prefer edge metadata if present.
    """
    if isinstance(edge, (list, tuple)) and len(edge) > 6:
        return [edge[6]]
    if gene_to_disorders is None:
        return []
    return [info[4] for info in gene_to_disorders.get(gene, [])]

# GRAND / CLUEreg helpers

_GRAND_RETRY = Retry(
    total=3,
    connect=3,
    read=3,
    backoff_factor=0.5,
    status_forcelist=(429, 500, 502, 503, 504),
    allowed_methods=("POST"),
    raise_on_status=False,
)

def _grand_session():
    s = requests.Session()
    adapter = HTTPAdapter(max_retries=_GRAND_RETRY)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    return s

def _grand_api_candidates(configured=None):
    base = (configured or os.getenv("GRAND_API", "https://grand.networkmedicine.org/api")).rstrip("/")
    fallbacks = [
        base,
        "https://grand.networkmedicine.org/api",
        "https://www.grand.networkmedicine.org/api",
    ]
    return list(dict.fromkeys(fallbacks))

def init_grand_api_context(grand_api=None, timeout=(3.05, 20)):
    return {
        "apis": _grand_api_candidates(grand_api),
        "timeout": timeout,
        "session": _grand_session(),
        "state": {
            "disabled": False,
            "error": "",
            "targeting_warning_reported": False,
            "cluereg_warning_reported": False,
        },
    }

def _grand_post(path, payload, context):
    state = context["state"]
    if state["disabled"]:
        raise RuntimeError(
            "GRAND API disabled for this run after previous failure.\n"
            + state["error"]
        )

    errors = []

    for base in context["apis"]:
        url = f"{base}{path}"
        try:
            r = context["session"].post(
                url,
                json=payload,
                timeout=context["timeout"],
            )
            r.raise_for_status()
            return r
        except requests.exceptions.RequestException as exc:
            errors.append(f"{url} :: {exc}")

    state["disabled"] = True
    state["error"] = (
        "failed to reach GRAND API; tried endpoints:\n"
        + "\n".join(errors)
    )
    raise RuntimeError(state["error"])

def local_targeting_scores(adj):
    # GRAND defines TF targeting as weighted outdegree.
    scores = {}
    for tf, edges in adj.items():
        score = 0.0
        for edge in edges:
            if len(edge) < 2:
                continue
            try:
                score += float(edge[1])
            except (TypeError, ValueError):
                continue
        scores[tf] = score
    return scores

def grand_targeting_scores(tf_list, context, tissue="brain"):
    """
    Query GRAND targeting scores via API.
    Returns dict {TF: targeting_score}.
    """
    if not tf_list:
        return {}

    payload = {
        "regulators": tf_list,
        "context": tissue,
        "species": "human",
    }

    r = _grand_post("/targeting", payload, context)
    data = r.json()
    scores = {}

    if isinstance(data, dict):
        data = data.get("results", [])

    for rec in data:
        tf = rec["regulator"]
        score = rec["targeting_score"]
        scores[tf] = score

    return scores

def differential_targeting_api(base_adj, disease_adj, context, tissue="brain"):
    # Compute differential targeting relative to healthy baseline.
    base_tfs = list(base_adj.keys())
    disease_tfs = list(disease_adj.keys())
    all_tfs = list(set(base_tfs) | set(disease_tfs))

    state = context["state"]
    score_source = "grand_api"

    try:
        base_scores = grand_targeting_scores(base_tfs, context, tissue=tissue)
        disease_scores = grand_targeting_scores(disease_tfs, context, tissue=tissue)
    except RuntimeError as exc:
        if not state["targeting_warning_reported"]:
            print("warning: GRAND targeting API unavailable; using local weighted outdegree fallback for all disorders")
            print(str(exc))
            state["targeting_warning_reported"] = True
        base_scores = local_targeting_scores(base_adj)
        disease_scores = local_targeting_scores(disease_adj)
        score_source = "local_fallback"

    diff_scores = []
    for tf in all_tfs:
        diff = disease_scores.get(tf, 0) - base_scores.get(tf, 0)
        diff_scores.append({
            "TF": tf,
            "diff_targeting": diff,
            "source": score_source,
        })

    diff_scores.sort(
        key=lambda x: abs(x["diff_targeting"]),
        reverse=True,
    )
    return diff_scores

def top_targeting_tfs(diff_scores, n=100):
    pos = [
        rec["TF"]
        for rec in diff_scores
        if rec["diff_targeting"] > 0
    ][:n]

    neg = [
        rec["TF"]
        for rec in diff_scores
        if rec["diff_targeting"] < 0
    ][:n]

    return pos, neg

def write_cluereg_file(tf_list, outfile):
    with open(outfile, "w") as f:
        for tf in tf_list:
            f.write(tf + "\n")

def cluereg_query(tf_pos, tf_neg, label, context, outfile):
    state = context["state"]
    if state["disabled"]:
        if not state["cluereg_warning_reported"]:
            print("warning: GRAND API unavailable; skipping CLUEreg API queries for all disorders")
            state["cluereg_warning_reported"] = True
        return

    if not tf_pos or not tf_neg:
        print(f"warning: skipping CLUEreg for {label}; need both non-empty positive and negative TF lists")
        return

    payload = {
        "up_regulators": tf_pos,
        "down_regulators": tf_neg,
    }

    try:
        r = _grand_post("/cluereg", payload, context)
    except RuntimeError as exc:
        print(f"warning: CLUEreg query failed for {label}")
        print(str(exc))
        return

    with open(outfile, "w") as f:
        f.write(r.text)

def build_cluereg_signatures_for_disorder(name, base_adj, disease_adj, context, outdir="results", tissue="brain", top_n=100):
    diff_scores = differential_targeting_api(
        base_adj,
        disease_adj,
        context,
        tissue=tissue,
    )
    pos, neg = top_targeting_tfs(diff_scores, n=top_n)

    write_cluereg_file(
        pos,
        f"{outdir}/cluereg_{name.lower()}_positive.txt",
    )
    write_cluereg_file(
        neg,
        f"{outdir}/cluereg_{name.lower()}_negative.txt",
    )
    cluereg_query(
        pos,
        neg,
        name.lower(),
        context,
        f"{outdir}/cluereg_results_{name.lower()}.json",
    )

    return {
        "diff_scores": diff_scores,
        "positive": pos,
        "negative": neg,
    }

# use below fns on deg_grn_both to return tf_scores
def regulatory_scores(deg_grn):

    """
    Calculate TF regulatory driver scores using PANDA edge weights
    and DEG log2FC values.

    driver_score = sum(weight * abs(log2fc))
    """

    tf_scores = []

    for tf, edges in deg_grn.items():

        score = 0
        targets = set()
        log2fcs = []

        for rec in edges:
            gene = rec[0]
            weight = rec[1]
            log2fc = rec[6]

            targets.add(gene)
            log2fcs.append(log2fc)

            score += weight * abs(log2fc)

        if not log2fcs:
            continue

        mean_fc = sum(log2fcs) / len(log2fcs)

        tf_scores.append({
            "TF": tf,
            "targets": len(targets),
            "edges": len(edges),
            "mean_log2fc": mean_fc,
            "driver_score": score})

    # sort a list of dictionaries by the value stored under "driver_score".
    tf_scores.sort(key=lambda x: x["driver_score"],reverse=True) 
    # TF with the largest driver score appears first.
    
    return tf_scores

# use below fns on deg_grn_tfsonly. 

def regulator_detection(grn, disorderlist=None): # outputs tf_results but loses information on sign of regulation.

    """
    Identify TFs whose high weight regulatory edges target genes with large expression shifts = regulator drivers in PANDA GRN network.

    Driver score per TF: sum(weight * abs(log2fc_target))
    = estimates how strongly a transcription factor's regulatory edges align with observed expression changes.
    aka total regulatory influence of a TF over the disease transcriptome
    
    Output: {'TF':..., 'targets':..., 'deg_targets':..., 'driver_score':..., 'mean_target_log2fc':...}
    """
    
    # Build gene: disorder metadata dictionary (only if needed)
    gene_to_disorders = None
    if disorderlist is not None:
        gene_to_disorders = {}
        for rec in disorderlist:
            gene_id = rec[0]
            info = tuple(rec[1:])
            gene_to_disorders.setdefault(gene_id, []).append(info)

    tf_results = []

    for tf, targets in grn.items():

        score = 0
        deg_targets = 0
        total_targets = 0
        log2fcs = []
        
        for edge in targets:
            gene = edge[0]
            weight = edge[1]
            total_targets += 1

            log2fc_vals = _edge_log2fc_vals(edge, gene, gene_to_disorders)
            if not log2fc_vals:
                continue

            deg_targets += 1

            for log2fc in log2fc_vals:
                score += weight * abs(log2fc)
                log2fcs.append(log2fc)
                
        if deg_targets == 0:
            continue

        mean_fc = sum(log2fcs) / len(log2fcs)
        
        normalized_score = score / total_targets 
        
        tf_results.append({
            "TF": tf,
            "targets": total_targets,
            "deg_targets": deg_targets,
            "driver_score": score, # = total regulatory influence
            "normalized_score": normalized_score, # = average influence per target
            "mean_target_log2fc": mean_fc})
    
    tf_results.sort(key=lambda x: x["driver_score"], reverse=True)

    return tf_results

def signed_regulator_detection(grn, disorderlist=None): # outputs tf_results. better than the last!

    """
    Detect regulatory drivers while preserving direction metadata.
    
    Identify TFs whose high weight regulatory edges target genes with large expression shifts = regulator drivers in PANDA GRN network.

    Driver score per TF: sum(weight * abs(log2fc_target))
    = estimates how strongly a transcription factor's regulatory edges align with observed expression changes.
    aka total regulatory influence of a TF over the disease transcriptome
    
    Output: {'TF':..., 'targets':..., 'deg_targets':..., 'driver_score':..., 'mean_target_log2fc':...}
    """

    gene_to_disorders = None
    if disorderlist is not None:
        gene_to_disorders = {}
        for rec in disorderlist:
            gene_id = rec[0]
            info = tuple(rec[1:])
            gene_to_disorders.setdefault(gene_id, []).append(info)

    tf_results = []

    for tf, targets in grn.items():

        driver_score = 0
        direction_score = 0

        deg_targets = 0
        total_targets = 0

        log2fcs = []
        up = 0
        down = 0

        for edge in targets:
            
            gene = edge[0]
            weight = edge[1]
            
            total_targets += 1

            log2fc_vals = _edge_log2fc_vals(edge, gene, gene_to_disorders)
            if not log2fc_vals:
                continue

            for log2fc in log2fc_vals:
                driver_score += weight * abs(log2fc)
                direction_score += weight * log2fc

                log2fcs.append(log2fc)

                deg_targets += 1

                if log2fc > 0:
                    up += 1
                elif log2fc < 0:
                    down += 1

        if deg_targets == 0:
            continue

        mean_fc = sum(log2fcs) / len(log2fcs)

        direction_ratio = direction_score / driver_score if driver_score != 0 else 0

        tf_results.append({
            "TF": tf,
            "targets": total_targets,
            "deg_targets": deg_targets,
            "driver_score": driver_score,
            "direction_score": direction_score,
            "direction_ratio": direction_ratio,
            "mean_target_log2fc": mean_fc,
            "up_targets": up,
            "down_targets": down
        })

    tf_results.sort(
        key=lambda x: x["driver_score"],
        reverse=True
    )

    return tf_results

def classify_drivers(tf_results):

    global_drivers = []
    focused_drivers = []

    for rec in tf_results:

        if rec["driver_score"] > 10 and rec["normalized_score"] < 0.05:
            global_drivers.append(rec)

        if rec["normalized_score"] > 0.08:
            focused_drivers.append(rec)

    return global_drivers, focused_drivers

def deg_enrichment(tf_results):

    enriched = []

    for rec in tf_results:

        fraction = rec["deg_targets"] / rec["targets"]

        if fraction > 0.15:
            enriched.append(rec)

    return enriched

def plot_regulators(tf_results, top_n = 10):

    """
    Plot top 10 TF regulatory drivers by driver_score
    """

    top = tf_results[:top_n]

    tfs = []
    scores = []

    for rec in top:
        tfs.append(rec["TF"])
        scores.append(rec["driver_score"])

    plt.figure(figsize=(10,6))

    plt.barh(tfs, scores)

    plt.xlabel("GRN Edge Weight x Differential Expression") # aka regulatory score
    plt.ylabel("Transcription Factor")

    plt.title("Top Regulatory Drivers")

    plt.gca().invert_yaxis()

    plt.tight_layout()
    plt.show()

def plot_deg_fraction(tf_results):

    """
    Scatter plot of TF targets vs DEG targets
    """

    targets = []
    deg_targets = []
    scores = []

    for rec in tf_results:

        targets.append(rec["targets"])
        deg_targets.append(rec["deg_targets"])
        scores.append(rec["driver_score"])

    plt.figure(figsize=(7,6))

    plt.scatter(targets, deg_targets)

    plt.xlabel("Total TF Targets")
    plt.ylabel("DEG Targets")

    plt.title("DEG enrichment across TF regulons")

    plt.tight_layout()
    plt.show()

def extract_modules(deggrn, tf_results, top_n = 10):

    modules = {}

    top_tfs = [i["TF"] for i in tf_results[:top_n]]

    for tf in top_tfs:

        targets = set()

        for i in deggrn[tf]:
            gene = i[0]
            targets.add(gene)

        modules[tf] = targets

    return modules

# use below fns for any adjlist of the structure

def edgeweight_summary(grn, degset):
    """
    Input
        adjlist : PANDA GRN adjacency list
                  {TF: [(target, weight), ...]}

        deg_set : set of differentially expressed genes

    Output
        dictionary summarizing edge weight statistics
    """

    tf_set = set(grn.keys()) 
    deg_set = _extract_deg_genes(degset)
    detf_set = tf_set & deg_set
    print("detfs",(list(detf_set)[:10]),list(deg_set)[:10])

    deg_weights = []
    nondeg_weights = []
    tf_means = []
    detf_means = []

    for tf, edges in grn.items():

        weights = []

        for edge in edges:
            
            gene = edge[0]
            weight = edge[1]
            
            weights.append(weight)

            if gene in deg_set:
                deg_weights.append(weight)
            else:
                nondeg_weights.append(weight)

        if weights:

            tf_mean = sum(weights) / len(weights)
            tf_means.append(tf_mean)

            if tf in detf_set:
                detf_means.append(tf_mean)

    # A large difference between DETF_mean_weight and TF_mean_weight MAY INDICATE some transcriptional driver behind a disorder's change
    
    summary = {
        "DEG_weight": sum(deg_weights)/len(deg_weights) if deg_weights else 0, # measures regulatory pressure acting on dysregulated genes. (TF to DEG)
        "nonDEG_weight": sum(nondeg_weights)/len(nondeg_weights) if nondeg_weights else 0, # baseline regulatory distribution (TF to background genes)
        "TF_weight": sum(tf_means)/len(tf_means) if tf_means else 0, # avg regulatory strength across all TF regulators.
        "DETF_weight": sum(detf_means)/len(detf_means) if detf_means else 0 # avg regulatory strength for DETFs
    }

    return summary

def edgeweight_summary_old(adjlist, *condition:None): # AVERAGES DOESN'T WORK
    """
    Input:
        Any adjacency list structured as:
            {TF: (gene1, regweight1), (gene2, regweight2)]}
    
    Output:
        printable dict containing average edge weight counts for:
              DETFs (TF nodes)
              DEGs (targets)
              DEXs (both tfs and gene targets)
            
    """
    tf_avg = [] # tf avgs
    g_avg = []
    degs_avg = [] # list of all DEG weights 
    detfs_avg = [] # no idea how to find this
    
    for tf, edges in adjlist.items():
        
        for edge in edges: # loops through each gene
            weight = edge[1]
            if len(edges) > 2:
                degs_avg.append(weight)
            else:
                g_avg.append(weight) # 
                
        tf_avg.append(sum(degs_avg) / len(degs_avg)) # adds average deg per tf
        
    summary = { 
        "DEGs_average":  sum(degs_avg) / len(degs_avg), # averages all degs
        "TFs_avg": sum(tf_avg) / len(tf_avg), # average TF
        "DETFs_avg": sum(),
        "G_avg": sum(g_avg) / len(g_avg)
        }
    
    return summary

# what about a function to look at differential TFs and the genes they regulate? 
# would require 

# DEG layer: aggregates all log2fc observations per gene. assigns gene sign by majority across observations
#TF layer: computes mean log2fc of DEG targets per TF. assigns TF “effect sign” based on its targets

def log2fc_summary(adjlist): # used chatgpt to fix previous function

    # gene  > list of log2fc values (prevents overwrite)
    gene_log2fc = {}

    # TF inferred effect sign (not DETF)
    tf_effect_sign = {}

    tf_avg = []   # per TF mean log2fc of its DEG targets
    deg_avg = []  # global DEG log2fc values

    for tf, edges in adjlist.items():

        tf_values = []   # per TF log2fc values
        tf_signs = []

        for edge in edges:

            # expect DEG enriched edge structure
            if len(edge) > 6:

                gene = edge[0]
                log2fc = edge[6]

                # store all observations per gene
                gene_log2fc.setdefault(gene, []).append(log2fc)

                deg_avg.append(log2fc)
                tf_values.append(log2fc)

                if log2fc > 0:
                    tf_signs.append("positive")
                elif log2fc < 0:
                    tf_signs.append("negative")
                else:
                    tf_signs.append("zero")

        # per TF average (only if it regulates ≥1 DEG)
        if tf_values:
            tf_mean = sum(tf_values) / len(tf_values)
            tf_avg.append(tf_mean)

        # infer TF effect direction from its DEG targets
        if tf_signs:
            pos = tf_signs.count("positive")
            neg = tf_signs.count("negative")

            if pos > neg:
                tf_effect_sign[tf] = "positive"
            elif neg > pos:
                tf_effect_sign[tf] = "negative"
            else:
                tf_effect_sign[tf] = "zero"

    # collapse gene level signs (majority across observations)
    deg_sign = {}
    for gene, vals in gene_log2fc.items():
        pos = sum(1 for v in vals if v > 0)
        neg = sum(1 for v in vals if v < 0)

        if pos > neg:
            deg_sign[gene] = "positive"
        elif neg > pos:
            deg_sign[gene] = "negative"
        else:
            deg_sign[gene] = "zero"

    def safe_mean(x):
        return sum(x) / len(x) if x else 0

    summary = {
        "DEGs_total": len(deg_sign),
        "DEGs_average": safe_mean(deg_avg),
        "DEGs_positive": sum(1 for s in deg_sign.values() if s == "positive"),
        "DEGs_negative": sum(1 for s in deg_sign.values() if s == "negative"),

        "TFs_total": len(tf_effect_sign),
        "TFs_avg": safe_mean(tf_avg),
        "TFs_positive": sum(1 for s in tf_effect_sign.values() if s == "positive"),
        "TFs_negative": sum(1 for s in tf_effect_sign.values() if s == "negative"),
    }

    return summary

# CSV-friendly summaries

def _safe_mean(vals):
    return sum(vals) / len(vals) if vals else 0

def _safe_median(vals):
    return float(np.median(vals)) if vals else 0

def deglist_stats(disorderlist):
    """
    Compute basic stats from a disorder list:
      [GENEID, DISORDER, STUDY, YEAR, TISSUE, LOG2FC, PVAL]
    """
    if not disorderlist:
        return {
            "records": 0,
            "unique_genes": 0,
            "mean_log2fc": 0,
            "median_log2fc": 0,
            "mean_abs_log2fc": 0,
            "pos": 0,
            "neg": 0,
            "zero": 0,
            "frac_pos": 0,
            "frac_neg": 0,
        }

    genes = [rec[0] for rec in disorderlist if rec]
    log2fcs = [rec[5] for rec in disorderlist if len(rec) > 5]

    pos = sum(1 for v in log2fcs if v > 0)
    neg = sum(1 for v in log2fcs if v < 0)
    zero = sum(1 for v in log2fcs if v == 0)
    total = len(log2fcs)

    return {
        "records": len(disorderlist),
        "unique_genes": len(set(genes)),
        "mean_log2fc": _safe_mean(log2fcs),
        "median_log2fc": _safe_median(log2fcs),
        "mean_abs_log2fc": _safe_mean([abs(v) for v in log2fcs]),
        "pos": pos,
        "neg": neg,
        "zero": zero,
        "frac_pos": pos / total if total else 0,
        "frac_neg": neg / total if total else 0,
    }

def adjlist_stats(adjlist):
    """
    Compute basic graph stats from an adjacency list.
    Supports both simple edges (gene, weight) and annotated edges.
    """
    if not adjlist:
        return {
            "tfs": 0,
            "genes": 0,
            "nodes": 0,
            "edge_rows": 0,
            "unique_edges": 0,
            "density": 0,
            "avg_out_degree": 0,
            "avg_in_degree": 0,
            "mean_weight": 0,
            "median_weight": 0,
            "mean_abs_log2fc": 0,
            "frac_up_log2fc": 0,
            "frac_down_log2fc": 0,
        }

    tfs = set(adjlist.keys())
    genes = set()
    weights = []
    log2fcs = []
    edge_rows = 0
    unique_edges = set()

    for tf, edges in adjlist.items():
        for edge in edges:
            if isinstance(edge, (list, tuple)):
                gene = edge[0]
                genes.add(gene)
                unique_edges.add((tf, gene))
                edge_rows += 1

                if len(edge) > 1:
                    try:
                        weights.append(float(edge[1]))
                    except (TypeError, ValueError):
                        pass

                if len(edge) > 6:
                    try:
                        log2fcs.append(float(edge[6]))
                    except (TypeError, ValueError):
                        pass
            else:
                genes.add(edge)
                unique_edges.add((tf, edge))
                edge_rows += 1

    nodes = tfs | genes
    n_tfs = len(tfs)
    n_genes = len(genes)
    unique_edge_count = len(unique_edges)

    density = unique_edge_count / (n_tfs * n_genes) if n_tfs and n_genes else 0
    avg_out = unique_edge_count / n_tfs if n_tfs else 0
    avg_in = unique_edge_count / n_genes if n_genes else 0

    pos = sum(1 for v in log2fcs if v > 0)
    neg = sum(1 for v in log2fcs if v < 0)
    total = len(log2fcs)

    return {
        "tfs": n_tfs,
        "genes": n_genes,
        "nodes": len(nodes),
        "edge_rows": edge_rows,
        "unique_edges": unique_edge_count,
        "density": density,
        "avg_out_degree": avg_out,
        "avg_in_degree": avg_in,
        "mean_weight": _safe_mean(weights),
        "median_weight": _safe_median(weights),
        "mean_abs_log2fc": _safe_mean([abs(v) for v in log2fcs]),
        "frac_up_log2fc": pos / total if total else 0,
        "frac_down_log2fc": neg / total if total else 0,
    }

# Bipartite GRN metrics (directed, weighted)

def _iter_weighted_edges(adjlist):
    """
    Yield (tf, gene, weight) triples from an adjacency list.
    Supports edge rows with extra metadata.
    """
    for tf, edges in adjlist.items():
        for edge in edges:
            if not isinstance(edge, (list, tuple)) or len(edge) < 2:
                continue
            gene = edge[0]
            try:
                weight = float(edge[1])
            except (TypeError, ValueError):
                continue
            yield tf, gene, weight

def grn_to_bipartite_digraph(adjlist, aggregate="max_abs"):
    """
    Build a directed bipartite graph from {TF: [(gene, weight, ...), ...]}.
    Duplicate TF->gene edges are aggregated to avoid duplicated disorder rows
    inflating graph measures.

    Important: TF and gene identifiers can overlap in GRNs (same symbol used
    as regulator and as target). To preserve strict bipartite structure for
    NetworkX bipartite algorithms, nodes are represented with typed IDs:
      ("TF", tf_id) and ("Gene", gene_id).
    """
    tf_nodes = set(adjlist.keys())
    gene_nodes = set()
    edge_weights = {}
    edge_counts = {}

    for tf, gene, weight in _iter_weighted_edges(adjlist):
        gene_nodes.add(gene)
        key = (tf, gene)
        edge_counts[key] = edge_counts.get(key, 0) + 1
        if key not in edge_weights:
            edge_weights[key] = weight
            continue
        if aggregate == "sum_abs":
            edge_weights[key] += abs(weight)
        elif aggregate == "mean":
            c = edge_counts[key]
            prev = edge_weights[key]
            edge_weights[key] = prev + (weight - prev) / c
        else:
            if abs(weight) > abs(edge_weights[key]):
                edge_weights[key] = weight

    G = nx.DiGraph()
    tf_node_ids = {("TF", tf) for tf in tf_nodes}
    gene_node_ids = {("Gene", gene) for gene in gene_nodes}

    for tf_id in tf_node_ids:
        G.add_node(tf_id, bipartite=0, kind="TF", label=tf_id[1])
    for gene_id in gene_node_ids:
        G.add_node(gene_id, bipartite=1, kind="Gene", label=gene_id[1])

    for (tf, gene), weight in edge_weights.items():
        abs_weight = abs(weight)
        if abs_weight <= 0:
            distance = float("inf")
        else:
            distance = 1.0 / abs_weight
        G.add_edge(
            ("TF", tf),
            ("Gene", gene),
            weight=weight,
            abs_weight=abs_weight,
            distance=distance,
        )
    return G, tf_node_ids, gene_node_ids

def _metric_partition(metric_map, tf_nodes, gene_nodes):
    tf_vals = [metric_map[n] for n in tf_nodes if n in metric_map]
    gene_vals = [metric_map[n] for n in gene_nodes if n in metric_map]
    return tf_vals, gene_vals

def _clean_numeric(values):
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return arr
    return arr[np.isfinite(arr)]

def weighted_node_degree(weighted_adjlist, node, mode="total"):
    """
    Weighted degree for directed bipartite GRNs.
      mode='out'   => outgoing abs-weight sum
      mode='in'    => incoming abs-weight sum
      mode='total' => in + out abs-weight sum
    """
    G, _, _ = grn_to_bipartite_digraph(weighted_adjlist)
    node_id = node
    if node_id not in G:
        tf_id = ("TF", node)
        gene_id = ("Gene", node)
        if tf_id in G and gene_id in G:
            node_id = None
        elif tf_id in G:
            node_id = tf_id
        elif gene_id in G:
            node_id = gene_id
        else:
            return 0.0

    def _node_val(n):
        if mode == "out":
            return float(G.out_degree(n, weight="abs_weight"))
        if mode == "in":
            return float(G.in_degree(n, weight="abs_weight"))
        return float(
            G.out_degree(n, weight="abs_weight")
            + G.in_degree(n, weight="abs_weight")
        )

    if node_id is None:
        # Identifier appears in both TF and Gene partitions; combine both roles.
        return _node_val(("TF", node)) + _node_val(("Gene", node))
    return _node_val(node_id)

def weighted_node_degree_stats(weighted_adjlist, mode="total"):
    G, _, _ = grn_to_bipartite_digraph(weighted_adjlist)
    vals = [
        weighted_node_degree(weighted_adjlist, node, mode=mode)
        for node in G.nodes
    ]
    return {
        "values": vals,
        "mean": _safe_mean(vals),
        "median": _safe_median(vals),
        "min": min(vals) if vals else 0.0,
        "max": max(vals) if vals else 0.0,
    }

def closeness_centralities(weighted_adj_list, mode="out"):
    """
    Weighted closeness using reciprocal edge-strength as distance.
    For directed graphs:
      mode='out' computes outward reachability centrality.
      mode='in' computes inward reachability centrality.
    """
    G, _, _ = grn_to_bipartite_digraph(weighted_adj_list)
    if G.number_of_nodes() == 0:
        return {}
    H = G.reverse(copy=False) if mode == "out" else G
    return nx.closeness_centrality(H, distance="distance", wf_improved=True)

def eccentricity_centralities(weighted_adj_list, mode="out"):
    """
    Reciprocal eccentricity (1/max shortest-path distance) with weighted edges.
    Unreachable-only nodes return 0.0.
    """
    G, _, _ = grn_to_bipartite_digraph(weighted_adj_list)
    if G.number_of_nodes() == 0:
        return {}
    H = G if mode == "out" else G.reverse(copy=False)
    ecc = {}
    for node in H.nodes:
        lengths = nx.single_source_dijkstra_path_length(
            H,
            node,
            weight="distance",
        )
        finite = [
            d for target, d in lengths.items()
            if target != node and np.isfinite(d)
        ]
        if not finite:
            ecc[node] = 0.0
        else:
            ecc[node] = 1.0 / max(finite)
    return ecc

def betweenness_centralities(weighted_adj_list, normalized=True):
    G, _, _ = grn_to_bipartite_digraph(weighted_adj_list)
    if G.number_of_nodes() == 0:
        return {}
    return nx.betweenness_centrality(
        G,
        normalized=normalized,
        weight="distance",
    )

def calculate_bipartite_clustering_coeff(weighted_adjlist):
    """
    Bipartite clustering for a TF-gene graph.
    Uses the bipartite overlap coefficient (networkx bipartite.clustering).
    """
    G, _, _ = grn_to_bipartite_digraph(weighted_adjlist)
    if G.number_of_nodes() == 0:
        return {}

    U = nx.Graph()
    U.add_nodes_from(G.nodes(data=True))
    U.add_edges_from((u, v) for u, v in G.edges())

    coeffs = nx_bipartite.clustering(U, mode="dot")
    for node in U.nodes():
        coeffs.setdefault(node, 0.0)
    return coeffs

def stat_test(group_a, group_b):
    """
    Mann-Whitney U test plus Cliff's delta effect size.
    """
    a = _clean_numeric(group_a)
    b = _clean_numeric(group_b)

    if a.size == 0 or b.size == 0:
        return {
            "U_statistic": float("nan"),
            "p_value": float("nan"),
            "cliffs_delta": float("nan"),
            "a_median": float("nan"),
            "b_median": float("nan"),
            "a_mean": float("nan"),
            "b_mean": float("nan"),
            "n_a": int(a.size),
            "n_b": int(b.size),
        }

    if mannwhitneyu is None:
        U = float("nan")
        p_val = float("nan")
    else:
        U, p_val = mannwhitneyu(a, b, alternative="two-sided")

    greater = 0
    less = 0
    for val in a:
        greater += np.sum(val > b)
        less += np.sum(val < b)
    cliff_delta = (greater - less) / (len(a) * len(b))

    return {
        "U_statistic": float(U) if np.isfinite(U) else U,
        "p_value": float(p_val) if np.isfinite(p_val) else p_val,
        "cliffs_delta": float(cliff_delta),
        "a_median": float(np.median(a)),
        "b_median": float(np.median(b)),
        "a_mean": float(np.mean(a)),
        "b_mean": float(np.mean(b)),
        "n_a": int(len(a)),
        "n_b": int(len(b)),
    }

def _probability_curve(values, bins=30):
    arr = _clean_numeric(values)
    if arr.size == 0:
        return np.array([]), np.array([])
    bin_count = min(bins, max(5, int(np.sqrt(arr.size) * 2)))
    counts, edges = np.histogram(arr, bins=bin_count)
    total = counts.sum()
    probs = counts / total if total else counts
    centers = (edges[:-1] + edges[1:]) / 2.0
    return centers, probs

def compare_dist(group_a, group_b, metric_name, outfile, label_a="TF", label_b="Gene"):
    """
    Plot two normalized metric distributions as probability curves.
    """
    x_a, y_a = _probability_curve(group_a)
    x_b, y_b = _probability_curve(group_b)

    fig, ax = plt.subplots(figsize=(8, 6))
    if x_a.size > 0:
        ax.plot(x_a, y_a, label=label_a, color="red")
    if x_b.size > 0:
        ax.plot(x_b, y_b, label=label_b, color="blue")
    ax.set_xlabel(metric_name)
    ax.set_ylabel("Probability")
    ax.set_title(f"{metric_name} Distribution Comparison")
    if ax.lines:
        ax.legend()
    else:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
    fig.savefig(outfile)
    plt.close(fig)

def bipartite_metric_distributions(weighted_adjlist):
    """
    Return run1-style metric distributions adapted for a directed,
    weighted bipartite GRN.
    """
    G, tf_nodes, gene_nodes = grn_to_bipartite_digraph(weighted_adjlist)
    if G.number_of_nodes() == 0:
        return {
            "Strength": {"TF": [], "Gene": []},
            "Clustering_Coefficient": {"TF": [], "Gene": []},
            "Closeness_Centrality": {"TF": [], "Gene": []},
            "Eccentricity_Centrality": {"TF": [], "Gene": []},
            "Betweenness": {"TF": [], "Gene": []},
        }

    strength_map = {}
    for node in G.nodes:
        strength_map[node] = float(
            G.in_degree(node, weight="abs_weight")
            + G.out_degree(node, weight="abs_weight")
        )

    clust_map = calculate_bipartite_clustering_coeff(weighted_adjlist)
    close_map = closeness_centralities(weighted_adjlist, mode="out")
    ecc_map = eccentricity_centralities(weighted_adjlist, mode="out")
    betw_map = betweenness_centralities(weighted_adjlist, normalized=True)

    tf_strength, gene_strength = _metric_partition(strength_map, tf_nodes, gene_nodes)
    tf_clust, gene_clust = _metric_partition(clust_map, tf_nodes, gene_nodes)
    tf_close, gene_close = _metric_partition(close_map, tf_nodes, gene_nodes)
    tf_ecc, gene_ecc = _metric_partition(ecc_map, tf_nodes, gene_nodes)
    tf_betw, gene_betw = _metric_partition(betw_map, tf_nodes, gene_nodes)

    return {
        "Strength": {"TF": tf_strength, "Gene": gene_strength},
        "Clustering_Coefficient": {"TF": tf_clust, "Gene": gene_clust},
        "Closeness_Centrality": {"TF": tf_close, "Gene": gene_close},
        "Eccentricity_Centrality": {"TF": tf_ecc, "Gene": gene_ecc},
        "Betweenness": {"TF": tf_betw, "Gene": gene_betw},
    }

def _plot_metric_panels(metric_data, disorder_name, boxplot_file, violinplot_file):
    order = [
        "Strength",
        "Clustering_Coefficient",
        "Closeness_Centrality",
        "Eccentricity_Centrality",
        "Betweenness",
    ]
    labels = {
        "Strength": "Strength",
        "Clustering_Coefficient": "Clust. Coeff.",
        "Closeness_Centrality": "Closeness",
        "Eccentricity_Centrality": "Eccentricity",
        "Betweenness": "Betweenness",
    }
    colors = ["lightcoral", "gold", "yellowgreen", "skyblue", "mediumpurple"]

    def _plot_ready(vals):
        arr = _clean_numeric(vals)
        if arr.size == 0:
            return [np.nan], False
        return arr.tolist(), True

    fig1, axs = plt.subplots(1, 5, figsize=(20, 4), layout="constrained")
    for i, metric in enumerate(order):
        tf_vals, tf_ok = _plot_ready(metric_data[metric]["TF"])
        gene_vals, gene_ok = _plot_ready(metric_data[metric]["Gene"])
        values = [tf_vals, gene_vals]
        if tf_ok or gene_ok:
            axs[i].boxplot(
                values,
                labels=["TF", "Gene"],
                patch_artist=True,
                boxprops={"facecolor": colors[i]},
            )
        else:
            axs[i].text(0.5, 0.5, "No data", ha="center", va="center", transform=axs[i].transAxes)
            axs[i].set_xticks([1, 2], ["TF", "Gene"])
        axs[i].set_title(labels[metric])
    fig1.suptitle(f"{disorder_name}: Boxplots of Bipartite Network Measures")
    fig1.savefig(boxplot_file)
    plt.close(fig1)

    fig2, axs = plt.subplots(1, 5, figsize=(20, 4), layout="constrained")
    for i, metric in enumerate(order):
        tf_vals, tf_ok = _plot_ready(metric_data[metric]["TF"])
        gene_vals, gene_ok = _plot_ready(metric_data[metric]["Gene"])
        values = [tf_vals, gene_vals]
        if tf_ok or gene_ok:
            axs[i].violinplot(values, showmeans=False, showmedians=True)
        else:
            axs[i].text(0.5, 0.5, "No data", ha="center", va="center", transform=axs[i].transAxes)
        axs[i].set_title(labels[metric])
        axs[i].set_xticks([1, 2], ["TF", "Gene"])
    fig2.suptitle(f"{disorder_name}: Violin Plots of Bipartite Network Measures")
    fig2.savefig(violinplot_file)
    plt.close(fig2)

def disorder_bipartite_report(adjlist, disorder_name, outdir="results"):
    """
    Compute run1-style statistics/visualizations for one disorder GRN.
    Returns a structured report with metric distributions, statistical tests,
    and output file paths.
    """
    metric_data = bipartite_metric_distributions(adjlist)
    os.makedirs(outdir, exist_ok=True)

    boxplot_file = os.path.join(outdir, f"{disorder_name.lower()}_bipartite_boxplots.png")
    violinplot_file = os.path.join(outdir, f"{disorder_name.lower()}_bipartite_violinplots.png")
    _plot_metric_panels(metric_data, disorder_name, boxplot_file, violinplot_file)

    stat_results = {}
    distribution_files = {}
    summary = {}

    for metric, groups in metric_data.items():
        tf_vals = groups["TF"]
        gene_vals = groups["Gene"]
        stat_results[metric] = stat_test(tf_vals, gene_vals)

        metric_slug = metric.lower()
        dist_file = os.path.join(outdir, f"{disorder_name.lower()}_{metric_slug}_distribution.png")
        compare_dist(tf_vals, gene_vals, metric, dist_file, label_a="TF", label_b="Gene")
        distribution_files[metric] = dist_file

        tf_arr = _clean_numeric(tf_vals)
        gene_arr = _clean_numeric(gene_vals)
        summary[f"{metric_slug}_tf_mean"] = float(np.mean(tf_arr)) if tf_arr.size else 0.0
        summary[f"{metric_slug}_gene_mean"] = float(np.mean(gene_arr)) if gene_arr.size else 0.0
        summary[f"{metric_slug}_tf_median"] = float(np.median(tf_arr)) if tf_arr.size else 0.0
        summary[f"{metric_slug}_gene_median"] = float(np.median(gene_arr)) if gene_arr.size else 0.0
        summary[f"{metric_slug}_p_value"] = stat_results[metric]["p_value"]
        summary[f"{metric_slug}_cliffs_delta"] = stat_results[metric]["cliffs_delta"]

    return {
        "metric_data": metric_data,
        "stats": stat_results,
        "summary": summary,
        "boxplot_file": boxplot_file,
        "violinplot_file": violinplot_file,
        "distribution_files": distribution_files,
    }

def disorder_bipartite_summary_row(name, report):
    """
    Flatten a disorder_bipartite_report into CSV-safe fields.
    """
    row = {"disorder": name}
    if not report:
        return row
    for key, value in report.get("summary", {}).items():
        row[f"bip_{key}"] = value
    row["bip_boxplot_file"] = report.get("boxplot_file", "")
    row["bip_violinplot_file"] = report.get("violinplot_file", "")
    return row

def disorder_bipartite_stat_rows(name, report):
    """
    Expanded per-metric table for one disorder's bipartite comparisons.
    """
    rows = []
    if not report:
        return rows
    for metric, stats in report.get("stats", {}).items():
        rows.append({
            "disorder": name,
            "metric": metric,
            "u_statistic": stats.get("U_statistic"),
            "p_value": stats.get("p_value"),
            "cliffs_delta": stats.get("cliffs_delta"),
            "tf_mean": stats.get("a_mean"),
            "gene_mean": stats.get("b_mean"),
            "tf_median": stats.get("a_median"),
            "gene_median": stats.get("b_median"),
            "n_tf": stats.get("n_a"),
            "n_gene": stats.get("n_b"),
        })
    return rows

def disorder_summary_row(
    name,
    disorderlist,
    detf_deggrn=None,
    tf_grn=None,
    deg_grn=None,
    reg_scores=None,
    tf_regulators=None,
    edgeweight_summary=None,
):
    """
    Build a single CSV-ready summary row for one disorder.
    """
    deg_stats = deglist_stats(disorderlist)

    row = {
        "disorder": name,
        "deg_records": deg_stats["records"],
        "deg_unique_genes": deg_stats["unique_genes"],
        "deg_mean_log2fc": deg_stats["mean_log2fc"],
        "deg_median_log2fc": deg_stats["median_log2fc"],
        "deg_mean_abs_log2fc": deg_stats["mean_abs_log2fc"],
        "deg_pos": deg_stats["pos"],
        "deg_neg": deg_stats["neg"],
        "deg_zero": deg_stats["zero"],
        "deg_frac_pos": deg_stats["frac_pos"],
        "deg_frac_neg": deg_stats["frac_neg"],
    }

    if detf_deggrn is not None:
        s = adjlist_stats(detf_deggrn)
        row.update({
            "detf_deggrn_tfs": s["tfs"],
            "detf_deggrn_targets": s["genes"],
            "detf_deggrn_edge_rows": s["edge_rows"],
            "detf_deggrn_unique_edges": s["unique_edges"],
            "detf_deggrn_density": s["density"],
            "detf_deggrn_mean_weight": s["mean_weight"],
            "detf_deggrn_median_weight": s["median_weight"],
            "detf_deggrn_mean_abs_log2fc": s["mean_abs_log2fc"],
        })

    if tf_grn is not None:
        s = adjlist_stats(tf_grn)
        row.update({
            "tf_grn_tfs": s["tfs"],
            "tf_grn_targets": s["genes"],
            "tf_grn_edge_rows": s["edge_rows"],
            "tf_grn_unique_edges": s["unique_edges"],
            "tf_grn_density": s["density"],
            "tf_grn_mean_weight": s["mean_weight"],
        })

    if deg_grn is not None:
        s = adjlist_stats(deg_grn)
        row.update({
            "deg_grn_tfs": s["tfs"],
            "deg_grn_targets": s["genes"],
            "deg_grn_edge_rows": s["edge_rows"],
            "deg_grn_unique_edges": s["unique_edges"],
            "deg_grn_density": s["density"],
            "deg_grn_mean_weight": s["mean_weight"],
        })

    if reg_scores is None and detf_deggrn is not None:
        reg_scores = regulatory_scores(detf_deggrn)
    if reg_scores:
        top = reg_scores[0]
        row["top_detf_driver"] = top.get("TF", "")
        row["top_detf_driver_score"] = top.get("driver_score", 0)
    else:
        row["top_detf_driver"] = ""
        row["top_detf_driver_score"] = 0

    if tf_regulators is None and tf_grn is not None:
        tf_regulators = regulator_detection(tf_grn, disorderlist)
    if tf_regulators:
        top = tf_regulators[0]
        row["top_tf_regulator"] = top.get("TF", "")
        row["top_tf_regulator_score"] = top.get("driver_score", 0)
    else:
        row["top_tf_regulator"] = ""
        row["top_tf_regulator_score"] = 0

    if edgeweight_summary:
        row["edgeweight_DEG_weight"] = edgeweight_summary.get("DEG_weight", 0)
        row["edgeweight_nonDEG_weight"] = edgeweight_summary.get("nonDEG_weight", 0)
        row["edgeweight_TF_weight"] = edgeweight_summary.get("TF_weight", 0)
        row["edgeweight_DETF_weight"] = edgeweight_summary.get("DETF_weight", 0)

    return row

def write_csv(rows, outfile, fieldnames=None):
    """
    Write list of dicts to CSV. Uses keys from first row if fieldnames not provided.
    """
    if not rows:
        return

    if fieldnames is None:
        fieldnames = list(rows[0].keys())

    with open(outfile, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

def pairwise_jaccard_rows(name_to_set):
    """
    Compute pairwise Jaccard similarity rows for a dict of sets.
    """
    names = sorted(name_to_set.keys())
    rows = []
    for i, a in enumerate(names):
        for b in names[i+1:]:
            aset = name_to_set[a]
            bset = name_to_set[b]
            inter = len(aset & bset)
            union = len(aset | bset)
            jaccard = inter / union if union else 0
            rows.append({
                "disorder_a": a,
                "disorder_b": b,
                "intersection": inter,
                "union": union,
                "jaccard": jaccard,
            })
    return rows

def log2fc_summary_old(adjlist): # only applies to degs = whatever has 7 values 
    #per tf (which does not have to be differentially regulated)
    
    """
    Input:
        DEG GRN adjacency list structured as :
            {TF: (DEG1, etc), (DEG2, etc)]}
            
        Figure out later if I should cut TFs down to DETF...
            {TF: (geneID, edgeweight, disorder, study, year, tissue, log2fc, pval)}

    Output:
        dict containing sign counts for:
              DEGs (targets)
              DETFs (TF nodes)
    """
    
    deg_sign = {}   # gene  > sign
    detf_sign = {}  # tf  > sign (relies on tf_signs within loop)
    tf_avg = []  # collects the avg reg weights of TFs that regulates a DEG 
    deg_avg = [] # collects the reg weights of DEGs
    # as calculated by averaging all its DEG target's log2fc   
    
    for tf, edges in adjlist.items(): # avgs might not work

        
        tf_signs = [] # collects all tf's pos or neg regulation per deg
        tf_weight = []
        for edge in edges: # loops through each gene
            
            if len(edge) > 2: # if differential gene, then length > 2
                weight = edge[6] # log2fc col
                gene = edge[0] # ensembl id col
                sign = ''
                if weight > 0.0:
                    sign = "positive"
                elif weight < 0.0:
                    sign = "negative"
                
                deg_sign[gene] = sign # key = DEG, value = sign
                tf_signs.append(sign) # adds sign to list of tf reg
                deg_avg.append(weight) 
                
        tf_avg += np.sum(deg_avg) / len(deg_avg) # averages all deg weights per tf
    
        if tf_signs: # tf differential values are obtained by comparing counts 
  
            pos = tf_signs.count("positive")
            neg = tf_signs.count("negative")

            if pos > neg:
                detf_sign[tf] = "positive"
            elif neg > pos:
                detf_sign[tf] = "negative"
            else:
                detf_sign[tf] = "zero"
                
    summary = { 
        "DEGs_total": len(deg_sign),
        "DEGs_average":  sum(deg_avg) / len(deg_avg), # averages all degs
        "DEGs_positive": sum(1 for s in deg_sign.values() if s == "positive"),
        "DEGs_negative": sum(1 for s in deg_sign.values() if s == "negative"),
        
        # Below are regular TFS THAT REGULATE A DEG (so maybe special tfs)
        "TFs_total": len(detf_sign),
        "TFs_avg": sum(tf_avg) / len(tf_avg), # average TF
        "TFs_positive": sum(1 for s in detf_sign.values() if s == "positive"),
        "TFs_negative": sum(1 for s in detf_sign.values() if s == "negative"),
        }

    return summary

def kmeans(adjlist): # struggles with clusters of different densities (aka evens out size)
    
    return

def calculate_clustering_coeff(adjlist): # input = adjlist like the one from A
    
    """
    Input:

    
    Output:  dict of (node:str, clustering_coeff:float)
    
    """
    C = {}
    clustering_coeff = 0
        
    for i in adjlist:
        
        # find total number of edges by adding all values of the input dictionary for denominator
        degree = len(adjlist[i])
        emax = degree * (degree-1)
        
        # check if denominator = 0 and equation is undefined
        if emax == 0:
            C[i] = 0.0 # because float
            continue # skip to next iteration of loop
        
        # numerator = check if neighbors of v are connected to each other
        
        edge = 0 # tracker
        # check every possible value of pair (u,w)
        for u in adjlist[i]:
            for w in adjlist[i]:
                #if w in adjlist[u]: # could also check if u in w's adj list for directed graphs         
                if u in adjlist and w in adjlist and w in adjlist[u]: # chat gpt corrected my original "if w in adjlist[u]"
                    edge += 1
        
        # formula 
        clustering_coeff = edge/emax # shouldn't it be edge/2 to correct for double counting?
        
        # C[i] is key
        C[i] = clustering_coeff
    
    return C

# I received help from CS drop in tutor Doran Penner for the below fn
def calculate_closeness(adjlist): #input adjlist from A, output dict of (node, closeness) pairs
    S = {} # dict to return
    closeness = 0.0
    #  closeness centrality = reciprocal of the sum of all nodes' shortest paths to v.
    # also see: len of adjlist   1 / sum of shortest paths to each nodes

    # initiate and populate distances dictionary
    
    nodelist = [] # all nodes excluding i 
    
    for i in adjlist:
        # 
        for n in adjlist: # make nodelist
            if n != i:
                nodelist.append(n)
    
        # length of shortest path from u to v = apply BFS          
        shortestpathdict = {} # shortest paths dictionary
        for k in nodelist:
            shortestpathdict[k] = float('inf')
            
        shortestpathdict[i] = 0 # set starting node to 0 
        
        # edited Anna's BFS code from lab2_sol.py to return shortest paths from s to all other nodes
        Q = [i] # initialize and populate queue of nodes to explore
        
        while len(Q) != 0:
            exploring = Q.pop(0)
            for neighbor in adjlist[exploring]: ## for each of exploring's neighbors...
                if shortestpathdict[neighbor] == float('inf'): # update if we haven't reached that neighbor yet
                    shortestpathdict[neighbor] = shortestpathdict[exploring] + 1
                    Q.append(neighbor)
        
        # hi: shortestpathdict bw non i and everything else
        ## calc closeness of i

        numerator = (len(adjlist))-1.0
        denominator = 0.0
        # with dictionary of distances, sum up all values from 'u's to v 
        denominator = sum(shortestpathdict.values())
    
        # apply formula
        closeness = numerator/denominator
        
        S[i] = closeness

    return S
