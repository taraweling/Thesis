import matplotlib.pyplot as plt
from pyvis.network import Network
import math
import numpy as np
import os
import re
import graph_utils as gu

# Shared process-level cache for ENSG -> display label lookups.
_GENE_LABEL_CACHE = {}
_GRN_POS_EDGE_COLOR = "#2e7d32"  # positive GRN weight
_GRN_NEG_EDGE_COLOR = "#7d3c98"  # negative GRN weight
_GRN_ZERO_EDGE_COLOR = "#656565"  # zero/unknown GRN weight
_LFC_POS_NODE_COLOR = "#4575b4"  # positive log2fc
_LFC_NEG_NODE_COLOR = "#d73027"  # negative log2fc
_LFC_NEUTRAL_NODE_COLOR = "#afaeae"


def _inject_disable_physics_after_stabilization(html_path):
    """
    Patch a PyVis HTML file so physics is disabled right after initial
    stabilization. This keeps a readable initial layout while improving
    interaction/load responsiveness afterward.
    """
    if not os.path.exists(html_path):
        return

    with open(html_path, "r", encoding="utf-8") as f:
        content = f.read()

    # Avoid duplicate injections if called multiple times.
    if "network.setOptions({ physics: { enabled: false } });" in content:
        return

    pattern = r'(network.once\("stabilizationIterationsDone", function\(\) \{)'
    replacement = (
        r'\1\n'
        r'                  network.setOptions({ physics: { enabled: false } });'
    )
    new_content, n = re.subn(pattern, replacement, content, count=1)
    if n == 0:
        return

    with open(html_path, "w", encoding="utf-8") as f:
        f.write(new_content)

def _write_empty_graph_html(outfile, message):
    """
    Overwrite outfile with a minimal HTML placeholder when no edges exist.
    This prevents stale prior visualizations from being mistaken for fresh output.
    """
    out_dir = os.path.dirname(outfile)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    html = (
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<title>Empty Graph</title></head><body>"
        f"<p>{message}</p>"
        "</body></html>"
    )
    with open(outfile, "w", encoding="utf-8") as f:
        f.write(html)


def _configure_network_runtime(network_obj, interactive, stabilization_iterations):
    """
    Apply runtime options for faster large-network rendering.
    """
    iterations = max(50, int(stabilization_iterations))
    configure_block = ""
    if interactive:
        configure_block = (
            '      "configure": {\n'
            '        "enabled": true,\n'
            '        "filter": ["physics"]\n'
            '      },\n'
        )
    options = f"""
    {{
{configure_block}
      "physics": {{
        "enabled": true,
        "stabilization": {{
          "enabled": true,
          "iterations": {iterations},
          "updateInterval": 50,
          "fit": true
        }}
      }},
      "interaction": {{
        "hideEdgesOnDrag": true,
        "hideEdgesOnZoom": true
      }}
    }}
    """
    network_obj.set_options(options)

def clear_gene_label_cache():
    """
    Clear cached gene labels (useful for tests/debugging).
    """
    _GENE_LABEL_CACHE.clear()

def _gene_display_label(gene_id):
    """
    Resolve a target gene node label for visualization.
    Keeps TF IDs unchanged; for ENSG targets, attempts symbol mapping and caches result.
    """
    if gene_id in _GENE_LABEL_CACHE:
        return _GENE_LABEL_CACHE[gene_id]

    if isinstance(gene_id, str) and gu.ENSEMBL_ID_RE.match(gene_id):
        mapped = gu.ensembl2gene(gene_id)
        label = mapped if mapped else gene_id
    else:
        label = gene_id

    _GENE_LABEL_CACHE[gene_id] = label
    return label

def _grn_edge_style(raw_weight, fallback_width=1.0):
    """
    Build edge width + color from GRN edge weight sign.
    Positive GRN weights are green, negative GRN weights are purple.
    """
    try:
        grn_weight = float(raw_weight)
    except (TypeError, ValueError):
        return fallback_width, _GRN_ZERO_EDGE_COLOR

    if grn_weight > 0:
        color = _GRN_POS_EDGE_COLOR
    elif grn_weight < 0:
        color = _GRN_NEG_EDGE_COLOR
    else:
        color = _GRN_ZERO_EDGE_COLOR

    return abs(grn_weight), color

def _lfc_node_color(log2fc_value):
    """
    Build node color from log2fc sign.
    Positive log2fc is blue, negative log2fc is red.
    """
    try:
        val = float(log2fc_value)
    except (TypeError, ValueError):
        return _LFC_NEUTRAL_NODE_COLOR

    if val > 0:
        return _LFC_POS_NODE_COLOR
    if val < 0:
        return _LFC_NEG_NODE_COLOR
    return _LFC_NEUTRAL_NODE_COLOR

def _merge_gene_node_metadata(node_entry, log2fc=None, pval=None):
    """
    Merge gene-level DEG metadata onto an existing node entry.
    Prefer the observation with the smallest p-value when available.
    """
    if not isinstance(node_entry, dict):
        return

    cur_pval = node_entry.get("pval")
    cur_log2fc = node_entry.get("log2fc")

    # Fill blanks first.
    if cur_log2fc is None and log2fc is not None:
        node_entry["log2fc"] = log2fc
    if cur_pval is None and pval is not None:
        node_entry["pval"] = pval
        if log2fc is not None:
            node_entry["log2fc"] = log2fc
        return

    # If both p-values are present, keep the stronger evidence (smaller p-value)
    # and carry its paired log2fc so color/size reflect the same observation.
    if cur_pval is not None and pval is not None:
        try:
            if float(pval) < float(cur_pval):
                node_entry["pval"] = pval
                if log2fc is not None:
                    node_entry["log2fc"] = log2fc
        except (TypeError, ValueError):
            pass

# All functions now assume edgelist rows look like:
# [TF, Gene, weight] OR
# [TF, Gene, weight, disorder, study, year, tissue, log2fc, pval]

# Visualize
## MUST FIGURE OUT A WAY TO VISUALIZE ONLY THE GENES THAT THE MOST TFS POINT TO 

#  returns graph of the most differentally expressed TFs 
def regulon_viz(adjlist, outfile ="regulon.html"): 
    return
    
def _top_tfs_adjlist(adjlist, top_tfs):
    if top_tfs is None:
        return adjlist

    tf_scores = {}
    for tf, edges in adjlist.items():
        targets = set()
        for edge in edges:
            if isinstance(edge, (list, tuple)):
                targets.add(edge[0])
            else:
                targets.add(edge)
        tf_scores[tf] = len(targets)

    top = sorted(tf_scores, key=tf_scores.get, reverse=True)[:top_tfs]
    return {tf: adjlist[tf] for tf in top}

def _top_degs_adjlist(adjlist, top_degs):
    if top_degs is None:
        return adjlist

    # Rank genes by |log2fc| desc, then pval asc (best evidence across edges)
    gene_best = {}
    for tf, edges in adjlist.items():
        for edge in edges:
            gene = edge[0] if isinstance(edge, (list, tuple)) else edge
            # adjlist edge format: (gene, weight, disorder, study, year, tissue, log2fc, pval)
            log2fc = edge[6] if isinstance(edge, (list, tuple)) and len(edge) > 6 else None
            pval = edge[7] if isinstance(edge, (list, tuple)) and len(edge) > 7 else None

            abs_log2fc = abs(log2fc) if log2fc is not None else 0.0
            pval_score = pval if pval is not None else float("inf")
            score = (-abs_log2fc, pval_score)

            if gene not in gene_best or score < gene_best[gene]:
                gene_best[gene] = score

    top = sorted(gene_best, key=gene_best.get)[:top_degs]
    keep = set(top)
    filtered = {}
    for tf, edges in adjlist.items():
        kept_edges = []
        for edge in edges:
            gene = edge[0] if isinstance(edge, (list, tuple)) else edge
            if gene in keep:
                kept_edges.append(edge)
        if kept_edges:
            filtered[tf] = kept_edges
    return filtered

def _top_degs_edgelist(edgelist, top_degs):
    if top_degs is None:
        return edgelist

    gene_best = {}
    for row in edgelist:
        gene = row[1]
        log2fc = row[7] if len(row) > 7 else None
        pval = row[8] if len(row) > 8 else None

        abs_log2fc = abs(log2fc) if log2fc is not None else 0.0
        pval_score = pval if pval is not None else float("inf")
        score = (-abs_log2fc, pval_score)

        if gene not in gene_best or score < gene_best[gene]:
            gene_best[gene] = score

    top = sorted(gene_best, key=gene_best.get)[:top_degs]
    keep = set(top)
    return [row for row in edgelist if row[1] in keep]

def _collapse_edgelist_tf_gene_pairs(edgelist):
    """
    Collapse repeated TF->gene rows in an edgelist before visualization.
    """
    adj = {}
    for row in edgelist:
        if len(row) < 3:
            continue
        tf = row[0]
        gene = row[1]
        edge = (gene, *row[2:])
        adj.setdefault(tf, []).append(edge)

    collapsed_adj = gu.aggregate_tf_gene_edges(adj, weight_agg="mean")

    collapsed = []
    for tf, edges in sorted(collapsed_adj.items(), key=lambda kv: str(kv[0])):
        for edge in edges:
            collapsed.append([tf, edge[0], *edge[1:]])
    return collapsed

def pyviz_deggrn(
    adjlist,
    outfile="grn.html",
    directed=True,
    top_tfs=None,
    top_degs=None,
    interactive=True,
    stabilization_iterations=300,
): # option to make not directed so it resembles a co expression net
    """
    Visualize a Gene Regulatory Network (GRN) stored as an adjacency list.

    Input formats supported:

    1) Simple:
        {TF: [(Gene, weight), ...]}

    2) Annotated:
        { TF: [(Gene, weight, disorder, study, year, tissue, log2fc, pval),
                ...]}

    Output:
        Interactive HTML visualization written to `outfile`.

    Visual encodings 
          Edge width   > |weight|
          Edge color   > GRN weight sign (green positive / purple negative)
          TF nodes     > triangle
          Gene nodes   > circle
          Node size    >  log10(pval)
          Node color   > log2fc (blue positive / red negative)
    Filtering:
          top_tfs  > keep TFs with most unique targets
          top_degs > keep genes with largest |log2fc|, then lowest pval
    """

    adjlist = gu.aggregate_tf_gene_edges(adjlist, weight_agg="mean")
    adjlist = _top_tfs_adjlist(adjlist, top_tfs)
    adjlist = _top_degs_adjlist(adjlist, top_degs)
    G = Network(directed=directed)

    # Collect node metadata

    node_info = {}

    for tf, targets in adjlist.items():
        node_info.setdefault(tf, {"type": "TF"})

        for edge in targets:
            gene = edge[0]
            weight = edge[1] if len(edge) > 1 else None

            # adjlist edge format: (gene, weight, disorder, study, year, tissue, log2fc, pval)
            log2fc = edge[6] if len(edge) > 6 else None
            pval = edge[7] if len(edge) > 7 else None

            gene_entry = node_info.setdefault(gene, {
                "type": "Gene",
                "log2fc": None,
                "pval": None,
            })
            _merge_gene_node_metadata(gene_entry, log2fc=log2fc, pval=pval)

    # Add nodes    

    for node, info in node_info.items():

        shape = "triangle" if info.get("type") == "TF" else "circle"
        label = node if info.get("type") == "TF" else _gene_display_label(node)

        # size from p value
        size = 20
        if info.get("pval") is not None and info["pval"] > 0:
            size = 10 + (-math.log10(info["pval"]) * 6)
            size = max(8, min(60, size))

        # color from log2fc
        color = _lfc_node_color(info.get("log2fc"))

        G.add_node(
            node,
            label=label,
            shape=shape,
            size=size,
            color=color,
        )

    # Add edges

    for tf, targets in adjlist.items():
        for edge in targets:
            gene = edge[0]
            raw_weight = edge[1] if len(edge) > 1 else 1.0
            weight, edge_color = _grn_edge_style(raw_weight, fallback_width=1.0)
            G.add_edge(tf, gene, value=weight, color=edge_color)

    # Write output
    if len(G.edges) == 0:
        print("Warning: pyviz_deggrn has no edges to render.")
        _write_empty_graph_html(outfile, "No edges to render for this graph.")
        return

    _configure_network_runtime(
        G,
        interactive=interactive,
        stabilization_iterations=stabilization_iterations,
    )
    G.write_html(outfile)
    return 

def viz_graph(
    edgelist,
    outfile,
    top_tfs=None,
    top_degs=None,
    interactive=True,
    stabilization_iterations=300,
): # constructed using anna ritz's course assignment for bio331
    """
    Visualize a directed graph from an enriched edgelist and write it to an HTML file.

    Input rows:
        [TF, Gene, weight, disorder, study, year, tissue, log2fc, pval]

    Visualization (if present):
          Edge width   > |weight|
          Edge color   > GRN weight sign (green positive / purple negative)
          TF nodes     > triangle
          Gene nodes   > circle
          Node size    >  log10(pval)
          Node color   > log2fc (blue positive / red negative)
    Filtering:
          top_tfs  > keep TFs with most targets
          top_degs > keep genes with largest |log2fc|, then lowest pval
    """

    edgelist = _collapse_edgelist_tf_gene_pairs(edgelist)

    if top_tfs is not None:
        tf_counts = {}
        for row in edgelist:
            tf = row[0]
            tf_counts[tf] = tf_counts.get(tf, 0) + 1
        top = set(sorted(tf_counts, key=tf_counts.get, reverse=True)[:top_tfs])
        edgelist = [row for row in edgelist if row[0] in top]
    edgelist = _top_degs_edgelist(edgelist, top_degs)

    G = Network(directed=True)

    # Collect nodes + metadata
    node_info = {}

    for row in edgelist:
        tf = row[0]
        gene = row[1]

        log2fc = row[7] if len(row) > 7 else None
        pval = row[8] if len(row) > 8 else None

        # TF node
        node_info.setdefault(tf, {"type": "TF"})

        # Gene node
        gene_entry = node_info.setdefault(gene, {
            "type": "Gene",
            "log2fc": None,
            "pval": None,
        })
        _merge_gene_node_metadata(gene_entry, log2fc=log2fc, pval=pval)

    # Add nodes with styling

    for n, info in node_info.items():

        # shape
        shape = "triangle" if info["type"] == "TF" else "circle"
        label = n if info["type"] == "TF" else _gene_display_label(n)

        # size from pval
        size = 20
        if info.get("pval") is not None and info["pval"] > 0:
            size = 10 + (-math.log10(info["pval"]) * 6)
            size = max(8, min(60, size))

        # color from log2fc
        color = _lfc_node_color(info.get("log2fc"))

        G.add_node(
            n,
            label=label,
            shape=shape,
            size=size,
            color=color,
        )

    # Add edge weight representations
    for row in edgelist:
        tf = row[0]
        gene = row[1]
        raw_weight = row[2] if len(row) > 2 else 1.0
        weight, edge_color = _grn_edge_style(raw_weight, fallback_width=1.0)
        G.add_edge(tf, gene, value=weight, color=edge_color)

    # Output
    if len(G.edges) == 0:
        print("Warning: viz_graph has no edges to render.")
        _write_empty_graph_html(outfile, "No edges to render for this graph.")
        return
    _configure_network_runtime(
        G,
        interactive=interactive,
        stabilization_iterations=stabilization_iterations,
    )
    G.write_html(outfile)

# Degree calculations
def get_degree(edgelist):
    """
    Returns a dictionary of (node, degree) pairs.

    Input rows:
        [TF, Gene, weight, ...metadata...]
    """

    degree = {}

    for row in edgelist:
        u = row[0]
        v = row[1]

        degree[u] = degree.get(u, 0) + 1
        degree[v] = degree.get(v, 0) + 1

    return degree

# IS THE BELOW CORRECT?
def calculate_degree(adjlist): # input = adj list
    #returns dict of (node:str, degree:int) pairs for every node in the graph 
    
    D = {} 
    degree = 0
    
    for i in adjlist:
        D[str(i)] = len(adjlist[i]) # explain this further

    return D

# Histogram utilities
def to_histogram(degree):
    """
    Converts {node: degree} into {degree: count}.
    """

    hist = {}

    for d in degree.values():
        hist[d] = hist.get(d, 0) + 1

    return hist

# Distribution plotting
def viz_distribution(hist1, hist2, outfile):
    """
    Plot two degree histograms and save to file.
    """

    x1, y1 = make_x_y(hist1)
    x2, y2 = make_x_y(hist2)

    fig, ax = plt.subplots()
    ax.plot(x1, y1, marker="o", label="G1")
    ax.plot(x2, y2, marker="s", label="G2")

    ax.set_xlabel("Degree")
    ax.set_ylabel("Number of Nodes")
    titlestr = 'Degree Distribution of' + str(hist1) + ' and ' + str(hist2)
    #ax.set_title("Degree Distribution of Graph1 and Graph2")
    ax.set_title(titlestr)
    ax.legend()

    plt.savefig(outfile)

# Helpers
def make_x_y(hist):
    x_list = []
    y_list = []
    for deg in sorted(hist.keys()):
        x_list.append(deg)
        y_list.append(hist[deg])
    return x_list, y_list


def _coerce_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def viz_bipartite_metric_stats(
    stat_rows,
    outfile,
    disorder_order=None,
    metric_order=None,
):
    """
    Visualize cross-disorder statistical comparisons for bipartite metrics.

    Input rows are expected from graph_algos.disorder_bipartite_stat_rows:
      {
        "disorder": ...,
        "metric": ...,
        "p_value": ...,
        "cliffs_delta": ...,
        ...
      }

    Produces a two-panel heatmap:
      - left: Cliff's delta (effect size, signed)
      - right: -log10(p-value) with significance stars
    """
    if not stat_rows:
        print("Warning: viz_bipartite_metric_stats has no rows to render.")
        return

    default_metric_order = [
        "Strength",
        "Clustering_Coefficient",
        "Closeness_Centrality",
        "Eccentricity_Centrality",
        "Betweenness",
    ]
    metric_order = metric_order or default_metric_order

    if disorder_order is None:
        disorder_order = sorted(
            {
                str(row.get("disorder", "")).strip()
                for row in stat_rows
                if row.get("disorder")
            }
        )

    disorder_order = [d for d in disorder_order if d]
    if not disorder_order:
        print("Warning: viz_bipartite_metric_stats found no disorders.")
        return

    metric_to_idx = {m: i for i, m in enumerate(metric_order)}
    disorder_to_idx = {d: i for i, d in enumerate(disorder_order)}

    cliffs_mat = np.full((len(metric_order), len(disorder_order)), np.nan)
    pval_mat = np.full((len(metric_order), len(disorder_order)), np.nan)
    sig_mat = np.zeros((len(metric_order), len(disorder_order)), dtype=bool)

    for row in stat_rows:
        metric = row.get("metric")
        disorder = row.get("disorder")
        if metric not in metric_to_idx or disorder not in disorder_to_idx:
            continue

        i = metric_to_idx[metric]
        j = disorder_to_idx[disorder]

        cliffs = _coerce_float(row.get("cliffs_delta"))
        pval = _coerce_float(row.get("p_value"))

        cliffs_mat[i, j] = cliffs
        pval_mat[i, j] = pval
        sig_mat[i, j] = np.isfinite(pval) and pval < 0.05

    # transform p-values for visibility while preserving NaN cells
    pscore_mat = np.full_like(pval_mat, np.nan, dtype=float)
    finite_mask = np.isfinite(pval_mat) & (pval_mat > 0)
    pscore_mat[finite_mask] = -np.log10(pval_mat[finite_mask])
    pscore_mat[np.isfinite(pval_mat) & (pval_mat == 0)] = 16.0

    fig, axs = plt.subplots(1, 2, figsize=(max(10, 1.4 * len(disorder_order) + 4), 7), layout="constrained")

    im1 = axs[0].imshow(cliffs_mat, cmap="coolwarm", vmin=-1, vmax=1, aspect="auto")
    axs[0].set_title("Cliff's Delta (TF vs Gene)")
    axs[0].set_xticks(range(len(disorder_order)))
    axs[0].set_xticklabels(disorder_order, rotation=45, ha="right")
    axs[0].set_yticks(range(len(metric_order)))
    axs[0].set_yticklabels(metric_order)
    plt.colorbar(im1, ax=axs[0], fraction=0.046, pad=0.04)

    max_pscore = np.nanmax(pscore_mat) if np.any(np.isfinite(pscore_mat)) else 1.0
    max_pscore = max(1.0, min(16.0, float(max_pscore)))
    im2 = axs[1].imshow(pscore_mat, cmap="YlOrRd", vmin=0, vmax=max_pscore, aspect="auto")
    axs[1].set_title("-log10(p-value)")
    axs[1].set_xticks(range(len(disorder_order)))
    axs[1].set_xticklabels(disorder_order, rotation=45, ha="right")
    axs[1].set_yticks(range(len(metric_order)))
    axs[1].set_yticklabels(metric_order)
    plt.colorbar(im2, ax=axs[1], fraction=0.046, pad=0.04)

    for i in range(len(metric_order)):
        for j in range(len(disorder_order)):
            if np.isfinite(cliffs_mat[i, j]):
                axs[0].text(j, i, f"{cliffs_mat[i, j]:.2f}", ha="center", va="center", fontsize=8)
            if np.isfinite(pscore_mat[i, j]):
                star = "*" if sig_mat[i, j] else ""
                axs[1].text(j, i, f"{pscore_mat[i, j]:.2f}{star}", ha="center", va="center", fontsize=8)

    fig.suptitle("Cross-Disorder Bipartite Metric Comparisons")
    fig.savefig(outfile, dpi=300)
    plt.close(fig)


def viz_strength_stat_summary(
    stat_rows,
    outfile,
    legend_outfile=None,
    disorder_order=None,
):
    """
    Create a Strength-only cross-disorder statistical summary figure.

    Input rows are expected from graph_algos.disorder_bipartite_stat_rows.
    Uses only rows where metric == "Strength".
    """
    if not stat_rows:
        print("Warning: viz_strength_stat_summary has no rows to render.")
        return

    strength_rows = [r for r in stat_rows if r.get("metric") == "Strength"]
    if not strength_rows:
        print("Warning: viz_strength_stat_summary found no Strength rows.")
        return

    rows_by_disorder = {r.get("disorder"): r for r in strength_rows if r.get("disorder")}
    if disorder_order is None:
        disorders = sorted(rows_by_disorder.keys())
    else:
        disorders = [d for d in disorder_order if d in rows_by_disorder]

    if not disorders:
        print("Warning: viz_strength_stat_summary found no disorders.")
        return

    cliffs = np.array([_coerce_float(rows_by_disorder[d].get("cliffs_delta")) for d in disorders], dtype=float)
    pvals = np.array([_coerce_float(rows_by_disorder[d].get("p_value")) for d in disorders], dtype=float)
    tf_means = np.array([_coerce_float(rows_by_disorder[d].get("tf_mean")) for d in disorders], dtype=float)
    gene_means = np.array([_coerce_float(rows_by_disorder[d].get("gene_mean")) for d in disorders], dtype=float)

    pscore = np.full_like(pvals, np.nan)
    valid_pos = np.isfinite(pvals) & (pvals > 0)
    pscore[valid_pos] = -np.log10(pvals[valid_pos])
    pscore[np.isfinite(pvals) & (pvals == 0)] = 16.0
    sig = np.isfinite(pvals) & (pvals < 0.05)

    colors = []
    for v in cliffs:
        if not np.isfinite(v):
            colors.append("#bdbdbd")
        elif v >= 0:
            colors.append("#c0392b")  # TF tends larger
        else:
            colors.append("#2471a3")  # Gene tends larger

    out_dir = os.path.dirname(outfile)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    x = np.arange(len(disorders))
    fig, axs = plt.subplots(
        2,
        1,
        figsize=(max(10, 1.2 * len(disorders) + 2), 9),
        layout="constrained",
        height_ratios=[2.0, 1.4],
    )

    # Panel A: effect size (Cliff's delta)
    axs[0].bar(x, cliffs, color=colors, edgecolor="black", linewidth=0.4)
    axs[0].axhline(0, color="black", linewidth=1)
    axs[0].set_ylim(-1.05, 1.05)
    axs[0].set_ylabel("Cliff's delta")
    axs[0].set_title("Panel A: Effect Size for Strength (TF vs Gene)")
    axs[0].set_xticks(x)
    axs[0].set_xticklabels(disorders, rotation=45, ha="right")
    for i, v in enumerate(cliffs):
        if np.isfinite(v):
            axs[0].text(i, v + (0.04 if v >= 0 else -0.06), f"{v:.2f}", ha="center", va="bottom" if v >= 0 else "top", fontsize=8)

    # Panel B: significance and group means
    max_pscore = np.nanmax(pscore) if np.any(np.isfinite(pscore)) else 1.0
    max_pscore = max(1.0, min(16.0, float(max_pscore)))
    axs[1].bar(x, pscore, color="#f4d03f", edgecolor="black", linewidth=0.4)
    axs[1].axhline(-np.log10(0.05), color="#7d3c98", linestyle="--", linewidth=1, label="p = 0.05")
    axs[1].set_ylim(0, max_pscore * 1.15)
    axs[1].set_ylabel("-log10(p-value)")
    axs[1].set_title("Panel B: Mann-Whitney U Significance for Strength")
    axs[1].set_xticks(x)
    axs[1].set_xticklabels(disorders, rotation=45, ha="right")
    for i, v in enumerate(pscore):
        if np.isfinite(v):
            star = "*" if sig[i] else ""
            axs[1].text(i, v + 0.03 * max(1.0, max_pscore), f"{v:.2f}{star}", ha="center", va="bottom", fontsize=8)
    axs[1].legend(loc="upper left")

    fig.suptitle("Cross-Disorder Strength Distribution Comparison")
    fig.savefig(outfile, dpi=300)
    plt.close(fig)

    legend_text = (
        "Figure legend: Cross-disorder statistical comparison of bipartite network Strength.\n"
        "Input: For each disorder, the DEG-filtered TF->gene GRN is converted to a directed bipartite graph. "
        "Node Strength is computed as weighted in-degree + weighted out-degree using absolute edge weights.\n"
        "Methods: TF-node and gene-node Strength distributions are compared within each disorder using a two-sided "
        "Mann-Whitney U test (non-parametric). Effect size is quantified by Cliff's delta, where positive values "
        "indicate TF strengths tend to exceed gene strengths, and negative values indicate the opposite.\n"
        "Panel A: Cliff's delta per disorder (magnitude and direction of effect).\n"
        "Panel B: -log10(p-value) from Mann-Whitney U; dashed line denotes p = 0.05; '*' marks p < 0.05.\n"
        "Additional context: Mean TF Strength and mean Gene Strength are also available in the stats table "
        "(columns tf_mean, gene_mean) used to construct this figure."
    )

    if legend_outfile:
        legend_dir = os.path.dirname(legend_outfile)
        if legend_dir:
            os.makedirs(legend_dir, exist_ok=True)
        with open(legend_outfile, "w") as f:
            f.write(legend_text + "\n")

    print(f"wrote strength stats figure: {outfile}")
    if legend_outfile:
        print(f"wrote strength stats legend: {legend_outfile}")
