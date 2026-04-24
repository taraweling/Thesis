import graph_utils as gu
import graph_viz as gv
import graph_algos as ga
import requests
import os
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

RUN_TESTS = False
TOP_DEGS = 10

GRAND_RETRY = Retry(
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
    adapter = HTTPAdapter(max_retries=GRAND_RETRY)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    return s


def _grand_api_candidates():
    configured = os.getenv("GRAND_API", "https://grand.networkmedicine.org/api").rstrip("/")
    fallbacks = [
        configured,
        "https://grand.networkmedicine.org/api",
        "https://www.grand.networkmedicine.org/api",
    ]
    # Deduplicate while preserving order.
    return list(dict.fromkeys(fallbacks))


# INSTRUCTIONS: RUN WHILE IN SRC FOLDER
def main():

    # obtain adjlist for each disorder (handle in main or helper fn?)
    ## what information do I want in my graph?
    """
    as a bipartite graph, node attribute should either be TF or gene.
    the TFs are DEGs and the target genes can be either DEGs or not differentially expressed genes in the GRAND regulatory network
    """

    # get the combined GRAND GRN (what check can I use?)
    brainother = gu.make_adjlist('data/Brain_Other.csv', 0.95)
    brainbg = gu.make_adjlist('data/Brain_Basal_Ganglia.csv', 0.95)

    # merge the two brain GRNs and test the merge behavior
    brains = gu.merge_adjlist(brainother, brainbg)
    if RUN_TESTS:
        test_merge_adjlist(brainother, brainbg, brains)
    grn = gu.ensemblify(brains)
    
    # produce text file of adjlist
    #with open("grand_brain_grn.txt", "w") as f: f.write("\n".join(f"{tf} [{gene}, {weight}]" for tf, edges in grn.items() for gene, weight in edges))
    gu.write_weighted_edgelist_txt(
        grn,
        "results/full_brain_grn_weighted_edgelist.txt",
    )
    full_ppi_path = os.path.join(gu.WALKER_INPUT_PPI_DIR, "full_brain_grn.ppi")
    gu.write_walker_ppi(grn, full_ppi_path)
    print("wrote Walker .ppi:", full_ppi_path)
    # run.py (right after ensemblify)
   
    ensg_tfs = [tf for tf in grn if gu.ENSEMBL_ID_RE.match(tf)]
    print("TFs mapped to ENSG:", len(ensg_tfs), "of", len(grn))

    # get adjlists per disorder by inputting the name of the disorder
    ## (options = AD, ADHD, BD, SZ, MDD, OCD) and the file location, outputting a list of lists
    disorders = ['AD', 'ADHD', 'ASD', 'BD', 'MDD', 'OCD', 'SZ']
    data_path = 'data/DEGDataStrictLFC.csv'
    degs = {d: gu.disorder_list(data_path, d) for d in disorders}

    # per-disorder outputs (the by-name suffix allows us to 
    detf_deggrn_by_name = {}
    detf_deggrn_viz_by_name = {}
    tf_grn_by_name = {}
    deg_grn_by_name = {}
    detf_deggrn_edgelist_by_name = {}
    reg_scores_by_name = {}
    tf_regulators_by_name = {}
    edgeweight_summary_by_name = {}
    log2fc_summary_by_name = {}
    bipartite_reports_by_name = {}
    summary_rows = []
    bipartite_stat_rows = []
    deg_gene_sets = {}
    detf_edge_sets = {}
    study_sets = {}
    tissue_sets = {}
    year_sets = {}

    for name, data in degs.items():
        if data is None:
            continue

        print("\n", "adjacency lists generated for", name)
        data = gu.ensemblifylist(data)
        degs[name] = data
        print("\n", "ensemblified DEGs for", name)

        degset = {row[0] for row in data}
        print(name, "deg first 10 keys:", list(degset)[:10])
        grn_tfs = set(grn.keys())
        grn_nodes = set(grn_tfs)
        for tf, edges in grn.items():
            for gene, _ in edges:
                grn_nodes.add(gene)
        print(name, "overlap (TFs):", len(grn_tfs & degset))
        print(name, "overlap (all GRN nodes):", len(grn_nodes & degset))

        # build DEG-filtered GRNs per disorder
        detf_deggrn = gu.de_grn_both(grn, data)  # DETFs and DEG targets
        tf_grn = gu.de_grn_tfsonly(grn, data)    # DETFs only
        deg_grn = gu.de_grn_degsonly(grn, data)  # DEG targets only

        # store by disorder name
        detf_deggrn_by_name[name] = detf_deggrn
        tf_grn_by_name[name] = tf_grn
        deg_grn_by_name[name] = deg_grn
        gu.write_weighted_edgelist_txt(
            detf_deggrn,
            f"results/{name.lower()}_detf_deggrn_weighted_edgelist.txt",
        )
        gu.write_weighted_edgelist_txt(
            deg_grn,
            f"results/{name.lower()}_deg_grn_weighted_edgelist.txt",
        )
        detf_ppi_path = os.path.join(
            gu.WALKER_INPUT_PPI_DIR, f"{name.lower()}_detf_deggrn.ppi"
        )
        deg_ppi_path = os.path.join(
            gu.WALKER_INPUT_PPI_DIR, f"{name.lower()}_deg_grn.ppi"
        )
        tf_ppi_path = os.path.join(
            gu.WALKER_INPUT_PPI_DIR, f"{name.lower()}_tf_grn.ppi"
        )
        gu.write_walker_ppi(detf_deggrn, detf_ppi_path)
        gu.write_walker_ppi(deg_grn, deg_ppi_path)
        gu.write_walker_ppi(tf_grn, tf_ppi_path)
        seed_txt_path = os.path.join(
            gu.WALKER_INPUT_SEED_DIR, f"{name.lower()}_seed.txt"
        )
        gu.save_adj_list_as_txt(deg_grn, seed_txt_path)
        print("wrote Walker .ppi:", detf_ppi_path)
        print("wrote Walker .ppi:", deg_ppi_path)
        print("wrote Walker .ppi:", tf_ppi_path)
        print("wrote Walker seed .txt:", seed_txt_path)

        # naming system for quick access (e.g., adhd_detf_deggrn)
        globals()[f"{name.lower()}_detf_deggrn"] = detf_deggrn
        globals()[f"{name.lower()}_tf_grn"] = tf_grn
        globals()[f"{name.lower()}_deg_grn"] = deg_grn

        # apply graph_algos to each disorder-specific GRN
        reg_scores_by_name[name] = ga.regulatory_scores(detf_deggrn)
        detf_deggrn_viz_by_name[name] = gu.aggregate_tf_gene_edges(detf_deggrn, weight_agg="mean")
        detf_deggrn_edgelist_by_name[name] = gu.adjlist2edgelist(detf_deggrn_viz_by_name[name])
        tf_regulators_by_name[name] = ga.regulator_detection(tf_grn, data)
        edgeweight_summary_by_name[name] = ga.edgeweight_summary(grn, data)
        log2fc_summary_by_name[name] = {
            "detf_deggrn": ga.log2fc_summary(detf_deggrn),
            "tf_grn": ga.log2fc_summary(tf_grn),
            "deg_grn": ga.log2fc_summary(deg_grn),
        }
        bipartite_reports_by_name[name] = ga.disorder_bipartite_report(
            detf_deggrn,
            name,
            outdir="results",
        )
        bipartite_stat_rows.extend(
            ga.disorder_bipartite_stat_rows(
                name,
                bipartite_reports_by_name[name],
            )
        )

        # collect per-disorder summary row for CSV
        disorder_row = ga.disorder_summary_row(
            name=name,
            disorderlist=data,
            detf_deggrn=detf_deggrn,
            tf_grn=tf_grn,
            deg_grn=deg_grn,
            reg_scores=reg_scores_by_name[name],
            tf_regulators=tf_regulators_by_name[name],
            edgeweight_summary=edgeweight_summary_by_name[name],
        )
        disorder_row.update(
            ga.disorder_bipartite_summary_row(
                name,
                bipartite_reports_by_name[name],
            )
        )
        summary_rows.append(disorder_row)

        # sets for pairwise Jaccard comparisons
        deg_gene_sets[name] = {rec[0] for rec in data}
        detf_edge_sets[name] = {
            (tf, edge[0])
            for tf, edges in detf_deggrn.items()
            for edge in edges
        }
        study_sets[name] = {rec[2] for rec in data if len(rec) > 2}
        tissue_sets[name] = {rec[4] for rec in data if len(rec) > 4}
        year_sets[name] = {rec[3] for rec in data if len(rec) > 3}

        def _edge_count(adj):
            return sum(len(edges) for edges in adj.values())

        def _top_regs(tf_scores, n=5):
            if not tf_scores:
                return "none"
            return ", ".join(
                f"{r['TF']}({r['driver_score']:.2f})"
                for r in tf_scores[:n]
            )

        print(f"\n{name} summary stats:")
        print(f"  detf_deggrn TFs={len(detf_deggrn)} edges={_edge_count(detf_deggrn)}")
        print(f"  tf_grn TFs={len(tf_grn)} edges={_edge_count(tf_grn)}")
        print(f"  deg_grn TFs={len(deg_grn)} edges={_edge_count(deg_grn)}")
        print(f"  top DETF regulators: {_top_regs(reg_scores_by_name[name])}")
        print(f"  top TF regulators: {_top_regs(tf_regulators_by_name[name])}")
        print(f"  edgeweight summary: {edgeweight_summary_by_name[name]}")
        print(f"  log2fc summary (detf_deggrn): {log2fc_summary_by_name[name]['detf_deggrn']}")
        print(f"  log2fc summary (tf_grn): {log2fc_summary_by_name[name]['tf_grn']}")
        print(f"  log2fc summary (deg_grn): {log2fc_summary_by_name[name]['deg_grn']}")
        print(f"  bipartite boxplot file: {bipartite_reports_by_name[name]['boxplot_file']}")
        print(f"  bipartite violin file: {bipartite_reports_by_name[name]['violinplot_file']}")

    # visualize graphs per disorder
    for name, detf_deggrn in detf_deggrn_by_name.items():
        detf_deggrn_viz = detf_deggrn_viz_by_name[name]
        detf_deggrnedgelist = detf_deggrn_edgelist_by_name[name]
        gv.viz_graph(
            detf_deggrnedgelist,
            f"results/deggrn_{name.lower()}.html",
            top_degs=TOP_DEGS,
        )
        gv.pyviz_deggrn(
            detf_deggrn_viz,
            outfile=f"results/pyviz_{name.lower()}.html",
            top_degs=TOP_DEGS,
        )

    # write CSV comparisons
    ga.write_csv(summary_rows, "results/deggrn_disorder_summary.csv")
    ga.write_csv(
        bipartite_stat_rows,
        "results/deggrn_bipartite_metric_stats.csv",
    )
    ga.write_csv(
        ga.pairwise_jaccard_rows(deg_gene_sets),
        "results/deggrn_jaccard_deg_genes.csv",)
    
    ga.write_csv(
        ga.pairwise_jaccard_rows(detf_edge_sets),
        "results/deggrn_jaccard_detf_edges.csv",)
    
    ga.write_csv(
        ga.pairwise_jaccard_rows(study_sets),
        "results/deggrn_jaccard_studies.csv",)
    
    ga.write_csv(
        ga.pairwise_jaccard_rows(tissue_sets),
        "results/deggrn_jaccard_tissues.csv",)
    
    ga.write_csv(
        ga.pairwise_jaccard_rows(year_sets),
        "results/deggrn_jaccard_years.csv",)
        
    """
    construct differential networks relative to a healthy baseline GRN
    compute TF differential targeting scores
    extract top 100 positively and negatively targeted TFs
    write CLUEreg-compatible lists
    """
    # Truong-style differential targeting using GRAND API
    # baseline healthy network = merged healthy brain GRN (Brain_Other + Basal_Ganglia + Cerebellum)
    print("\nquerying GRAND API for CLUEreg-compatible TF signatures")
    GRAND_APIS = _grand_api_candidates()
    GRAND_TIMEOUT = (3.05, 20)
    GRAND_SESSION = _grand_session()
    GRAND_API_STATE = {
        "disabled": False,
        "error": "",
        "targeting_warning_reported": False,
        "cluereg_warning_reported": False,
    }

    def _grand_post(path, payload):
        if GRAND_API_STATE["disabled"]:
            raise RuntimeError(
                "GRAND API disabled for this run after previous failure.\n"
                + GRAND_API_STATE["error"]
            )

        errors = []

        for base in GRAND_APIS:
            url = f"{base}{path}"
            try:
                r = GRAND_SESSION.post(
                    url,
                    json=payload,
                    timeout=GRAND_TIMEOUT,
                )
                r.raise_for_status()
                return r
            except requests.exceptions.RequestException as exc:
                errors.append(f"{url} :: {exc}")

        GRAND_API_STATE["disabled"] = True
        GRAND_API_STATE["error"] = (
            "failed to reach GRAND API; tried endpoints:\n"
            + "\n".join(errors)
        )
        raise RuntimeError(GRAND_API_STATE["error"])

    def _local_targeting_scores(adj):
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

    def _grand_targeting_scores(tf_list, tissue="brain"):
        """
        query GRAND targeting scores via API
        returns dict {TF: targeting_score}
        """

        if not tf_list:
            return {}

        payload = {
            "regulators": tf_list,
            "context": tissue,
            "species": "human"}

        r = _grand_post("/targeting", payload)

        data = r.json()

        scores = {}

        if isinstance(data, dict):
            # Some API responses are wrapped in an object.
            data = data.get("results", [])

        for rec in data:

            tf = rec["regulator"]
            score = rec["targeting_score"]

            scores[tf] = score

        return scores

    def _differential_targeting_api(base_adj, disease_adj):
        #compute differential targeting relative to healthy baseline

        base_tfs = list(base_adj.keys())
        disease_tfs = list(disease_adj.keys())

        all_tfs = list(set(base_tfs) | set(disease_tfs))

        score_source = "grand_api"
        try:
            base_scores = _grand_targeting_scores(base_tfs)
            disease_scores = _grand_targeting_scores(disease_tfs)
        except RuntimeError as exc:
            if not GRAND_API_STATE["targeting_warning_reported"]:
                print("warning: GRAND targeting API unavailable; using local weighted outdegree fallback for all disorders")
                print(str(exc))
                GRAND_API_STATE["targeting_warning_reported"] = True
            base_scores = _local_targeting_scores(base_adj)
            disease_scores = _local_targeting_scores(disease_adj)
            score_source = "local_fallback"

        diff_scores = []

        for tf in all_tfs:

            diff = disease_scores.get(tf, 0) - base_scores.get(tf, 0)

            diff_scores.append({
                "TF": tf,
                "diff_targeting": diff,
                "source": score_source})

        diff_scores.sort(
            key=lambda x: abs(x["diff_targeting"]),
            reverse=True)

        return diff_scores


    def _top100(diff_scores):

        pos = [
            rec["TF"]
            for rec in diff_scores
            if rec["diff_targeting"] > 0][:100]

        neg = [
            rec["TF"]
            for rec in diff_scores
            if rec["diff_targeting"] < 0][:100]

        return pos, neg


    def _write_cluereg_file(tf_list, outfile):

        with open(outfile, "w") as f:

            for tf in tf_list:
                f.write(tf + "\n")


    def _cluereg_query(tf_pos, tf_neg, label):
        if GRAND_API_STATE["disabled"]:
            if not GRAND_API_STATE["cluereg_warning_reported"]:
                print("warning: GRAND API unavailable; skipping CLUEreg API queries for all disorders")
                GRAND_API_STATE["cluereg_warning_reported"] = True
            return

        if not tf_pos or not tf_neg:
            print(f"warning: skipping CLUEreg for {label}; need both non-empty positive and negative TF lists")
            return

        payload = {
            "up_regulators": tf_pos,
            "down_regulators": tf_neg}

        try:
            r = _grand_post("/cluereg", payload)
        except RuntimeError as exc:
            print(f"warning: CLUEreg query failed for {label}")
            print(str(exc))
            return

        with open(f"results/cluereg_results_{label}.json", "w") as f:

            f.write(r.text)


    clue_lists = {}
    diff_scores_by_name = {}


    for name in disorders:

        if name not in detf_deggrn_by_name:
            continue

        print("\nprocessing disorder:", name)

        disease_net = detf_deggrn_by_name[name]

        diff_scores = _differential_targeting_api(
            grn,
            disease_net)

        diff_scores_by_name[name] = diff_scores

        pos100, neg100 = _top100(diff_scores)

        clue_lists[name] = {
            "positive": pos100,
            "negative": neg100}

        _write_cluereg_file(
            pos100,
            f"results/cluereg_{name.lower()}_positive.txt")

        _write_cluereg_file(
            neg100,
            f"results/cluereg_{name.lower()}_negative.txt")

        _cluereg_query(
            pos100,
            neg100,
            name.lower())

        print(
            name,
            "CLUEreg TFs:",
            len(pos100),
            "positive,",
            len(neg100),
            "negative")

    _, _, ffl_rows, fbl_rows = gu.detect_regulatory_loops(
        degs,
        detf_deggrn_by_name,
    )

    ga.write_csv(
        ffl_rows,
        "results/deggrn_feedforward_loops.csv",
    )
    ga.write_csv(
        fbl_rows,
        "results/deggrn_feedback_loops.csv",
    )


def test_merge_adjlist(brainother, brainbg, brains):
    def _edge_set(adj):
        return {(tf, gene) for tf, edges in adj.items() for gene, _ in edges}

    def _edge_weights(adj):
        return {(tf, gene): w for tf, edges in adj.items() for gene, w in edges}

    edges_other = _edge_set(brainother)
    edges_bg = _edge_set(brainbg)
    edges_merged = _edge_set(brains)

    # 1) merged contains all edges from both inputs
    assert edges_other.issubset(edges_merged)
    assert edges_bg.issubset(edges_merged)
    assert edges_merged == edges_other | edges_bg

    # 2) overlap is on edges (TF, Gene), not just TF keys
    edge_overlap = edges_other & edges_bg
    print("merge successful", "brain grn edge overlap:", len(edge_overlap))

    # 3) merged weights are averaged for overlapping edges (spot check one)
    if edge_overlap:
        sample = next(iter(edge_overlap))
        w1 = _edge_weights(brainother)[sample]
        w2 = _edge_weights(brainbg)[sample]
        wm = _edge_weights(brains)[sample]
        assert abs(wm - ((w1 + w2) / 2)) < 1e-9

if __name__ == '__main__':
    main()
