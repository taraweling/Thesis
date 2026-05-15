import graph_utils as gu
import graph_viz as gv
import graph_algos as ga
import requests
import os
import re
import csv
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

RUN_TESTS = False
TOP_DEGS = 10
TOP_TFS = 10
VIZ_INTERACTIVE = True
VIZ_STABILIZATION_ITERATIONS = 300 # controls how long PyVis physics runs during the initial layout stabilization of each graph.
GRN_ZSCORE_THRESHOLD = 1.96
GRAND_TARGETING_MODE = os.getenv("GRAND_TARGETING_MODE", "local").strip().lower()
CLUE_SIGNATURE_MODE = os.getenv("CLUE_SIGNATURE_MODE", "gene").strip().lower()
if CLUE_SIGNATURE_MODE not in {"tf", "gene"}:
    CLUE_SIGNATURE_MODE = "tf"
CLUE_STRICT_SIGNED = os.getenv("CLUE_STRICT_SIGNED", "0").strip().lower() not in {
    "0",
    "false",
    "no",
    "off",
}
try:
    CLUE_TOP_N = max(1, int(os.getenv("CLUE_TOP_N", "100")))
except ValueError:
    CLUE_TOP_N = 100
TOPOLOGY_TOP_N = 10
TOPOLOGY_METRIC_SPECS = [
    ("strength", "strength"),
    ("clustering", "clustering_coefficient"),
    ("closeness", "closeness_centrality"),
    ("eccentricity", "eccentricity_centrality"),
    ("betweenness", "betweenness"),
]

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
        "https://grand.networkmedicine.org/api/v1",
        "https://www.grand.networkmedicine.org/api/v1",
    ]
    # Deduplicate while preserving order.
    return list(dict.fromkeys(fallbacks))


def _grand_analysis_candidates():
    configured = os.getenv(
        "GRAND_ANALYSIS_URL",
        "https://grand.networkmedicine.org/analysis/",
    ).rstrip("/")
    fallbacks = [
        configured,
        "https://grand.networkmedicine.org/analysis",
        "https://www.grand.networkmedicine.org/analysis",
    ]
    # Normalize with trailing slash and deduplicate while preserving order.
    normalized = []
    for url in fallbacks:
        if not url.endswith("/"):
            url = url + "/"
        normalized.append(url)
    return list(dict.fromkeys(normalized))

_ENSEMBL_STABLE_RE = re.compile(r"^ENSG(\d{11})$")


def _repair_probable_broken_ensembl_ids(records, valid_target_ids):
    """
    Repair likely mangled ENSG IDs seen in some DEG sources.

    Heuristic:
      - candidate must look like ENSG + 11 digits
      - candidate currently does not match a known GRN target
      - candidate starts with ENSG00001/ENSG00002 (suspicious in this dataset)
      - drop one trailing digit from the significant numeric part (likely version
        suffix leak), zero-pad back to 11, and accept only if it exists in GRN
        targets
    """
    repaired = []
    changed = 0
    recovered = 0

    for rec in records:
        if not rec:
            continue

        gene_id = rec[0]
        new_gene_id = gene_id

        if gene_id not in valid_target_ids:
            m = _ENSEMBL_STABLE_RE.match(gene_id)
            if m and (gene_id.startswith("ENSG00001") or gene_id.startswith("ENSG00002")):
                sig = m.group(1).lstrip("0") or "0"
                if len(sig) > 1:
                    candidate = "ENSG" + sig[:-1].zfill(11)
                    if candidate in valid_target_ids:
                        new_gene_id = candidate

        if new_gene_id != gene_id:
            changed += 1
            if new_gene_id in valid_target_ids:
                recovered += 1

        repaired.append([new_gene_id, *rec[1:]])

    return repaired, changed, recovered


def _repair_disorder_records(records, valid_target_ids):
    """
    Apply conservative ENSG repair heuristic to one disorder's DEG records.
    Returns repaired records and diagnostic stats.
    """
    records = records or []
    raw_degset = {rec[0] for rec in records if rec}
    raw_target_overlap = len(valid_target_ids & raw_degset)
    low_overlap_threshold = max(10, int(0.01 * max(1, len(raw_degset))))

    stats = {
        "repair_applied": False,
        "changed_records": 0,
        "recovered_records": 0,
        "raw_target_overlap": raw_target_overlap,
        "repaired_target_overlap": raw_target_overlap,
    }

    repaired_records = records
    if raw_target_overlap < low_overlap_threshold:
        candidate_records, changed_count, recovered_count = _repair_probable_broken_ensembl_ids(
            records,
            valid_target_ids,
        )
        repaired_degset = {rec[0] for rec in candidate_records if rec}
        repaired_overlap = len(valid_target_ids & repaired_degset)
        if repaired_overlap > raw_target_overlap:
            repaired_records = candidate_records
            stats.update(
                {
                    "repair_applied": True,
                    "changed_records": changed_count,
                    "recovered_records": recovered_count,
                    "repaired_target_overlap": repaired_overlap,
                }
            )

    return repaired_records, stats


def _write_repaired_degdata_csv(
    input_csv_path,
    output_csv_path,
    disorders,
    valid_target_ids,
):
    """
    Write one all-disorder CSV with repaired ENSG IDs, preserving DEGData.csv
    column order/shape exactly.
    """
    with open(input_csv_path, newline="", encoding="utf-8-sig") as infile:
        reader = csv.reader(infile)
        all_rows = list(reader)

    if not all_rows:
        print("warning: DEG input CSV is empty; skipped repaired DEG CSV export")
        return

    header = all_rows[0]
    data_rows = all_rows[1:]
    try:
        disorder_idx = header.index("DISORDER")
        gene_idx = header.index("GENEID")
    except ValueError:
        print(
            "warning: could not find DISORDER/GENEID columns; "
            "skipped repaired DEG CSV export"
        )
        return

    disorder_set = set(disorders)
    row_indices_by_disorder = {d: [] for d in disorders}
    records_by_disorder = {d: [] for d in disorders}

    for i, row in enumerate(data_rows):
        if disorder_idx >= len(row):
            continue
        disorder_name = row[disorder_idx].strip()
        if disorder_name not in disorder_set:
            continue
        gene_id = row[gene_idx].strip() if gene_idx < len(row) else ""
        row_indices_by_disorder[disorder_name].append(i)
        records_by_disorder[disorder_name].append([gene_id])

    total_changed = 0
    total_recovered = 0
    applied_disorders = 0

    for disorder_name in disorders:
        disorder_records = records_by_disorder.get(disorder_name, [])
        repaired_records, stats = _repair_disorder_records(disorder_records, valid_target_ids)
        if stats["repair_applied"]:
            applied_disorders += 1
            total_changed += stats["changed_records"]
            total_recovered += stats["recovered_records"]
            print(
                f"{disorder_name} DEG CSV repair: changed {stats['changed_records']} "
                f"records, recovered {stats['recovered_records']} GRN target matches "
                f"({stats['raw_target_overlap']} -> {stats['repaired_target_overlap']})"
            )

        for row_i, repaired_rec in zip(row_indices_by_disorder[disorder_name], repaired_records):
            if not repaired_rec:
                continue
            if gene_idx >= len(data_rows[row_i]):
                data_rows[row_i].extend([""] * (gene_idx - len(data_rows[row_i]) + 1))
            data_rows[row_i][gene_idx] = repaired_rec[0]

    with open(output_csv_path, "w", newline="", encoding="utf-8") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(header)
        writer.writerows(data_rows)

    print(
        "wrote repaired DEG csv:",
        output_csv_path,
        f"(repair applied in {applied_disorders} disorders;",
        f"records changed={total_changed}; recovered={total_recovered})",
    )


# INSTRUCTIONS: RUN WHILE IN SRC FOLDER
def main():

    # obtain adjlist for each disorder (handle in main or helper fn?)
    ## what information do I want in my graph?
    """
    as a bipartite graph, node attribute should either be TF or gene.
    the TFs are DEGs and the target genes can be either DEGs or not differentially expressed genes in the GRAND regulatory network
    """

    # get the combined GRAND GRN (what check can I use?)
    brainother = gu.make_adjlist('data/Brain_Other.csv', GRN_ZSCORE_THRESHOLD)
    brainbg = gu.make_adjlist('data/Brain_Basal_Ganglia.csv', GRN_ZSCORE_THRESHOLD)
    braincerebellum = gu.make_adjlist('data/Brain_Cerebellum.csv', GRN_ZSCORE_THRESHOLD)

    # merge the brain GRNs and test the merge behavior
    brains = gu.merge_adjlist(brainother, brainbg, braincerebellum)
    if RUN_TESTS:
        test_merge_adjlist(brainother, brainbg, braincerebellum, brains)
    grn = gu.ensemblify_targets(brains)
    
    # produce text file of adjlist
    #with open("grand_brain_grn.txt", "w") as f: f.write("\n".join(f"{tf} [{gene}, {weight}]" for tf, edges in grn.items() for gene, weight in edges))
    gu.write_weighted_edgelist_txt(
        grn,
        "results/full_brain_grn_weighted_edgelist.txt",
    )
    # PPI export disabled per request.
    # full_ppi_path = os.path.join(gu.WALKER_INPUT_PPI_DIR, "full_brain_grn.ppi")
    # gu.write_walker_ppi(grn, full_ppi_path)
    # print("wrote Walker .ppi:", full_ppi_path)
    # TF/source nodes remain symbols; target genes are mapped to ENSG where possible.
    ensg_tfs = [tf for tf in grn if gu.ENSEMBL_ID_RE.match(tf)]
    target_nodes = {gene for edges in grn.values() for gene, _ in edges}
    ensg_targets = [gene for gene in target_nodes if gu.ENSEMBL_ID_RE.match(gene)]
    # print("TFs mapped to ENSG:", len(ensg_tfs), "of", len(grn))
    print("targets mapped to ENSG:", len(ensg_targets), "of", len(target_nodes))
    grn_target_nodes = set(target_nodes)

    # get adjlists per disorder by inputting the name of the disorder
    ## (options = AD, ADHD, BD, SZ, MDD, OCD) and the file location, outputting a list of lists
    disorders = ['AD', 'ADHD', 'ASD', 'BD', 'MDD', 'OCD', 'SZ']
    data_path = 'data/DEGData.csv'
    repaired_deg_csv_path = "results/DEGData_repaired_ENSG_all_disorders.csv"
    _write_repaired_degdata_csv(
        data_path,
        repaired_deg_csv_path,
        disorders,
        grn_target_nodes,
    )
    degs = {d: gu.disorder_list(data_path, d) for d in disorders}

    # per-disorder outputs
    tf_grn_by_name = {}
    deg_grn_by_name = {}
    detf_deggrn_by_name = {}
    network_variants_by_name = {}
    reg_scores_by_name = {}
    tf_regulators_by_name = {}
    edgeweight_summary_by_name = {}
    log2fc_summary_by_name = {}
    variant_labels = ("deggrn", "detfdeggrn", "detfgrn")
    network_variants_by_name = {
        disorder: {label: {} for label in variant_labels}
        for disorder in disorders
    }
    bipartite_reports_by_variant = {label: {} for label in variant_labels}
    bipartite_stat_rows_by_variant = {label: [] for label in variant_labels}
    summary_rows = []
    deg_gene_sets = {}
    tf_edge_sets = {}
    study_sets = {}
    tissue_sets = {}
    year_sets = {}

    for name, data in degs.items():
        if data is None:
            print(f"warning: no DEG rows for {name}; generating empty graph outputs")
            continue

        print("\n", "adjacency lists generated for", name)
        raw_data = data  # keep original DEG Gene ID symbols for TF overlap checks
        raw_data, repair_stats = _repair_disorder_records(raw_data, grn_target_nodes)
        raw_degset = {row[0] for row in raw_data if row}
        if repair_stats["repair_applied"]:
            print(
                f"{name} DEG ID repair: changed {repair_stats['changed_records']} records, "
                f"recovered {repair_stats['recovered_records']} GRN target matches "
                f"({repair_stats['raw_target_overlap']} -> {repair_stats['repaired_target_overlap']})"
            )

        data = gu.ensemblifylist(raw_data)  # target gene processes remain ENSG-based
        degs[name] = data
        print("\n", "ensemblified DEGs for", name)

        degset = {row[0] for row in data}
        print(name, "deg first 10 keys:", list(degset)[:10])
        grn_tfs = set(grn.keys())
        grn_nodes = set()
        for tf, edges in grn.items():
            for gene, _ in edges:
                grn_nodes.add(gene)
        print(name, "overlap (TFs by DEG Gene ID):", len(grn_tfs & raw_degset))
        print(name, "overlap (target ENSG nodes):", len(grn_nodes & degset))

        # build DEG-filtered GRNs per disorder using mixed IDs:
        # TF/source overlap by raw DEG Gene ID symbols, targets by ENSG IDs.
        tf_grn = gu.de_grn_detfgrn(
            grn,
            data,
            tf_degset=raw_degset,
            strict_tf_filter=False,
        )
        deg_grn = gu.de_grn_degsonly(grn, data, tf_degset=raw_degset)
        
        # Strict DETF-DEG network: both TF and target must be DE.
        detf_deggrn = gu.de_grn_tfsonly(
            grn,
            data,
            tf_degset=raw_degset,
            strict_tf_filter=True,
        )

        # store by disorder name
        tf_grn_by_name[name] = tf_grn
        deg_grn_by_name[name] = deg_grn
        detf_deggrn_by_name[name] = detf_deggrn
        network_variants = {
            "deggrn": deg_grn,
            "detfdeggrn": detf_deggrn,
            "detfgrn": tf_grn,
        }
        network_variants_by_name[name] = network_variants
        gu.write_weighted_edgelist_txt(
            tf_grn,
            f"results/{name.lower()}_tf_grn_weighted_edgelist.txt",
        )
        gu.write_weighted_edgelist_txt(
            deg_grn,
            f"results/{name.lower()}_deg_grn_weighted_edgelist.txt",
        )
        # Per-disorder PPI/seed exports disabled per request.
        # deg_ppi_path = os.path.join(
        #     gu.WALKER_INPUT_PPI_DIR, f"{name.lower()}_deg_grn.ppi"
        # )
        # tf_ppi_path = os.path.join(
        #     gu.WALKER_INPUT_PPI_DIR, f"{name.lower()}_tf_grn.ppi"
        # )
        # gu.write_walker_ppi(deg_grn, deg_ppi_path)
        # gu.write_walker_ppi(tf_grn, tf_ppi_path)
        # seed_txt_path = os.path.join(
        #     gu.WALKER_INPUT_SEED_DIR, f"{name.lower()}_seed.txt"
        # )
        # gu.save_adj_list_as_txt(deg_grn, seed_txt_path)
        # print("wrote Walker .ppi:", deg_ppi_path)
        # print("wrote Walker .ppi:", tf_ppi_path)
        # print("wrote Walker seed .txt:", seed_txt_path)

        # naming system for quick access
        globals()[f"{name.lower()}_tf_grn"] = tf_grn
        globals()[f"{name.lower()}_deg_grn"] = deg_grn

        # apply graph_algos to each disorder-specific GRN
        reg_scores_by_name[name] = ga.regulatory_scores(detf_deggrn)
        tf_regulators_by_name[name] = ga.regulator_detection(tf_grn, data)
        edgeweight_summary_by_name[name] = ga.edgeweight_summary(
            grn,
            data,
            tf_degset=raw_degset,
        )
        log2fc_summary_by_name[name] = {
            "detf_deggrn": ga.log2fc_summary(detf_deggrn),
            "tf_grn": ga.log2fc_summary(tf_grn),
            "deg_grn": ga.log2fc_summary(deg_grn),
        }
        for variant, variant_adj in network_variants.items():
            report = ga.disorder_bipartite_report(
                variant_adj,
                f"{variant}_{name}",
                outdir="results",
            )
            bipartite_reports_by_variant[variant][name] = report
            bipartite_stat_rows_by_variant[variant].extend(
                ga.disorder_bipartite_stat_rows(name, report)
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
                bipartite_reports_by_variant["detfgrn"].get(name, {}),
            )
        )
        summary_rows.append(disorder_row)

        # sets for pairwise Jaccard comparisons
        deg_gene_sets[name] = {rec[0] for rec in data}
        tf_edge_sets[name] = {
            (tf, edge[0])
            for tf, edges in tf_grn.items()
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
        print(f"  top TF regulators (weighted by DEG shifts): {_top_regs(reg_scores_by_name[name])}")
        print(f"  top TF regulators: {_top_regs(tf_regulators_by_name[name])}")
        print(f"  edgeweight summary: {edgeweight_summary_by_name[name]}")
        print(f"  log2fc summary (detf_deggrn): {log2fc_summary_by_name[name]['detf_deggrn']}")
        print(f"  log2fc summary (tf_grn): {log2fc_summary_by_name[name]['tf_grn']}")
        print(f"  log2fc summary (deg_grn): {log2fc_summary_by_name[name]['deg_grn']}")
        for variant in variant_labels:
            report = bipartite_reports_by_variant[variant].get(name, {})
            print(f"  {variant} bipartite boxplot file: {report.get('boxplot_file', '')}")
            print(f"  {variant} bipartite violin file: {report.get('violinplot_file', '')}")

    # Re-render bipartite box/violin panels with shared y-limits per metric
    # so values are directly comparable across disorders, within each variant.
    for variant in variant_labels:
        shared_y_limits = ga.metric_panel_y_limits(
            [
                report.get("metric_data")
                for report in bipartite_reports_by_variant[variant].values()
                if report
            ]
        )
        bipartite_stat_rows_by_variant[variant] = []
        for name, network_variants in network_variants_by_name.items():
            prior_report = bipartite_reports_by_variant[variant].get(name, {})
            report = ga.disorder_bipartite_report(
                network_variants.get(variant, {}),
                f"{variant}_{name}",
                outdir="results",
                y_limits=shared_y_limits,
                metric_data=prior_report.get("metric_data"),
            )
            bipartite_reports_by_variant[variant][name] = report
            bipartite_stat_rows_by_variant[variant].extend(
                ga.disorder_bipartite_stat_rows(name, report)
            )

    # Keep healthy panels on the same metric axes as disorder panels.
    disorder_metric_data = [
        report.get("metric_data")
        for variant in variant_labels
        for report in bipartite_reports_by_variant[variant].values()
        if report
    ]
    shared_disorder_y_limits = ga.metric_panel_y_limits(disorder_metric_data)
    healthy_bipartite_report = ga.disorder_bipartite_report(
        grn,
        "healthy_brain_grn",
        outdir="results",
        y_limits=shared_disorder_y_limits,
    )
    healthy_bipartite_stat_rows = ga.disorder_bipartite_stat_rows(
        "Healthy_Brain_GRN",
        healthy_bipartite_report,
    )

    # visualize graphs per disorder and network variant
    for name, network_variants in network_variants_by_name.items():
        for label, net_adj in network_variants.items():
            net_viz = gu.aggregate_tf_gene_edges(net_adj, weight_agg="mean")
            net_edgelist = gu.adjlist2edgelist(net_viz)

            gv.viz_graph(
                net_edgelist,
                f"results/{label}_{name.lower()}.html",
                top_tfs=TOP_TFS,
                top_degs=TOP_DEGS,
                interactive=VIZ_INTERACTIVE,
                stabilization_iterations=VIZ_STABILIZATION_ITERATIONS,
            )
            gv.pyviz_deggrn(
                net_viz,
                outfile=f"results/pyviz_{label}_{name.lower()}.html",
                top_tfs=TOP_TFS,
                top_degs=TOP_DEGS,
                interactive=VIZ_INTERACTIVE,
                stabilization_iterations=VIZ_STABILIZATION_ITERATIONS,
            )
    expected_html = []
    for name in disorders:
        for label in variant_labels:
            expected_html.append(f"results/{label}_{name.lower()}.html")
            expected_html.append(f"results/pyviz_{label}_{name.lower()}.html")
    missing_html = [path for path in expected_html if not os.path.exists(path)]
    if missing_html:
        print("warning: missing graph html outputs:")
        for path in missing_html:
            print(" ", path)
    else:
        print(
            "confirmed graph html outputs:",
            f"{len(expected_html)} files across {len(disorders)} disorders and {len(variant_labels)} graph types",
        )

    # write CSV comparisons
    ga.write_csv(summary_rows, "results/deggrn_disorder_summary.csv")
    combined_bipartite_stat_rows = []
    for variant in variant_labels:
        variant_stat_rows = bipartite_stat_rows_by_variant[variant]
        ga.write_csv(
            variant_stat_rows,
            f"results/{variant}_bipartite_metric_stats.csv",
        )
        combined_bipartite_stat_rows.extend(
            [{"graph_type": variant, **row} for row in variant_stat_rows]
        )
        gv.viz_bipartite_metric_stats(
            variant_stat_rows,
            f"results/{variant}_bipartite_metric_stats_heatmap.png",
            disorder_order=disorders,
        )
        gv.viz_strength_stat_summary(
            variant_stat_rows,
            f"results/{variant}_strength_stats_summary.png",
            legend_outfile=f"results/{variant}_strength_stats_figure_legend.txt",
            disorder_order=disorders,
        )
        print(
            "wrote bipartite metric comparison plot:",
            f"results/{variant}_bipartite_metric_stats_heatmap.png",
        )

    ga.write_csv(
        healthy_bipartite_stat_rows,
        "results/healthy_brain_grn_bipartite_metric_stats.csv",
    )
    combined_bipartite_stat_rows.extend(
        [{"graph_type": "healthy_brain_grn", **row} for row in healthy_bipartite_stat_rows]
    )
    if combined_bipartite_stat_rows:
        ga.write_csv(
            combined_bipartite_stat_rows,
            "results/combined_bipartite_metric_stats.csv",
            fieldnames=[
                "graph_type",
                "disorder",
                "metric",
                "u_statistic",
                "p_value",
                "cliffs_delta",
                "tf_mean",
                "gene_mean",
                "tf_median",
                "gene_median",
                "n_tf",
                "n_gene",
            ],
        )
        print(
            "wrote combined bipartite stats csv:",
            "results/combined_bipartite_metric_stats.csv",
        )
    print(
        "wrote healthy GRN bipartite boxplot:",
        healthy_bipartite_report.get("boxplot_file", ""),
    )
    ga.write_csv(
        ga.pairwise_jaccard_rows(deg_gene_sets),
        "results/deggrn_jaccard_deg_genes.csv",)
    
    ga.write_csv(
        ga.pairwise_jaccard_rows(tf_edge_sets),
        "results/deggrn_jaccard_tf_edges.csv",)
    
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
    GRAND_ANALYSIS_URLS = _grand_analysis_candidates()
    GRAND_TIMEOUT = (3.05, 20)
    GRAND_SESSION = _grand_session()
    GRAND_API_STATE = {
        "disabled": False,
        "error": "",
        "targeting_warning_reported": False,
        "cluereg_warning_reported": False,
        "analysis_warning_reported": False,
        "targeting_mode_reported": False,
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

    def _local_gene_targeting_scores(adj):
        # GRAND defines gene targeting as weighted indegree.
        scores = {}
        for _, edges in adj.items():
            for edge in edges:
                if len(edge) < 2:
                    continue
                gene = edge[0]
                try:
                    weight = float(edge[1])
                except (TypeError, ValueError):
                    continue
                scores[gene] = scores.get(gene, 0.0) + weight
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

        use_api = GRAND_TARGETING_MODE in {"api", "grand_api"}
        score_source = "local_fallback"

        if use_api:
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
        else:
            if not GRAND_API_STATE["targeting_mode_reported"]:
                print(
                    "note: GRAND_TARGETING_MODE=local; "
                    "using local weighted outdegree targeting for all disorders"
                )
                GRAND_API_STATE["targeting_mode_reported"] = True
            base_scores = _local_targeting_scores(base_adj)
            disease_scores = _local_targeting_scores(disease_adj)

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

    def _differential_gene_targeting_local(base_adj, disease_adj):
        base_scores = _local_gene_targeting_scores(base_adj)
        disease_scores = _local_gene_targeting_scores(disease_adj)
        all_genes = set(base_scores) | set(disease_scores)

        diff_scores = []
        for gene in all_genes:
            diff = disease_scores.get(gene, 0.0) - base_scores.get(gene, 0.0)
            diff_scores.append({
                "Gene": gene,
                "diff_targeting": diff,
                "source": "local_fallback",
            })

        diff_scores.sort(key=lambda x: abs(x["diff_targeting"]), reverse=True)
        return diff_scores

    def _mean_log2fc_by_gene(disorder_records):
        sums = {}
        counts = {}
        for rec in disorder_records:
            if len(rec) < 6:
                continue
            gene = rec[0]
            try:
                log2fc = float(rec[5])
            except (TypeError, ValueError):
                continue
            sums[gene] = sums.get(gene, 0.0) + log2fc
            counts[gene] = counts.get(gene, 0) + 1
        return {gene: sums[gene] / counts[gene] for gene in sums if counts.get(gene)}


    def _top_n_signature(
        diff_scores,
        key_name="TF",
        n=100,
        strict_signed=True,
        allowed_values=None,
    ):
        """
        Build CLUEreg high-targeted / low-targeted lists from differential
        targeting scores.

        Preferred behavior:
          - up list: strictly positive differential targeting
          - down list: strictly negative differential targeting

        Optional fallback behavior (strict_signed=False):
          - up list: top-ranked TFs by differential targeting (most positive / least negative)
          - down list: bottom-ranked TFs by differential targeting (most negative / least positive)
        """
        if not diff_scores:
            return [], []

        if allowed_values is not None:
            allowed_values = set(allowed_values)
            diff_scores = [
                rec for rec in diff_scores
                if rec.get(key_name) in allowed_values
            ]
            if not diff_scores:
                return [], []

        by_desc = sorted(diff_scores, key=lambda x: x["diff_targeting"], reverse=True)
        by_asc = sorted(diff_scores, key=lambda x: x["diff_targeting"])

        # Non-strict mode: always rank by differential targeting and return
        # top-N / bottom-N regardless of sign balance.
        if not strict_signed:
            pos = [rec[key_name] for rec in by_desc][:n]
            neg = [rec[key_name] for rec in by_asc][:n]
            return pos, neg

        pos = [rec[key_name] for rec in by_desc if rec["diff_targeting"] > 0][:n]
        neg = [rec[key_name] for rec in by_asc if rec["diff_targeting"] < 0][:n]

        if strict_signed or (pos and neg):
            return pos, neg

        # Rank-based fallback so CLUEreg can still receive two directional lists.
        if not pos:
            used = set(neg)
            pos = [rec[key_name] for rec in by_desc if rec[key_name] not in used][:n]
        if not neg:
            used = set(pos)
            neg = [rec[key_name] for rec in by_asc if rec[key_name] not in used][:n]

        return pos, neg


    def _write_cluereg_file(tf_list, outfile):

        with open(outfile, "w") as f:

            for tf in tf_list:
                f.write(tf + "\n")


    def _cluereg_signature_filename(label, entity_label, direction, n):
        """
        Build a clear, descriptive signature filename.
        Example:
          cluereg_ad_top100_most_positive_differentially_targeted_genes.txt
        """
        return (
            "results/"
            f"cluereg_{label}_top{n}_{direction}_"
            f"differentially_targeted_{entity_label}.txt"
        )


    def _tf_topology_metric_map(adjlist, metric_name="betweenness", deg_gene_candidates=None):
        G, tf_nodes, gene_nodes = ga.grn_to_bipartite_digraph(adjlist)
        if not tf_nodes:
            return {}

        if metric_name == "strength":
            metric_map = {
                node: float(
                    G.in_degree(node, weight="abs_weight")
                    + G.out_degree(node, weight="abs_weight")
                )
                for node in G.nodes
            }
        elif metric_name == "closeness":
            metric_map = ga.closeness_centralities(adjlist, mode="out")
        elif metric_name == "eccentricity":
            metric_map = ga.eccentricity_centralities(adjlist, mode="out")
        elif metric_name == "clustering":
            metric_map = ga.calculate_bipartite_clustering_coeff(adjlist)
        elif metric_name == "betweenness":
            # TF betweenness over shortest paths BETWEEN DEG target-gene pairs.
            # Use an undirected bipartite graph so gene<->TF<->gene paths exist.
            U = G.to_undirected()
            if deg_gene_candidates is None:
                deg_gene_nodes = sorted(gene_nodes)
            else:
                deg_gene_nodes = sorted(
                    {
                        ("Gene", gene)
                        for gene in deg_gene_candidates
                        if ("Gene", gene) in U
                    }
                )

            metric_map = dict.fromkeys(U.nodes, 0.0)
            if len(deg_gene_nodes) >= 2:
                raw_map = ga.nx.betweenness_centrality_subset(
                    U,
                    sources=deg_gene_nodes,
                    targets=deg_gene_nodes,
                    normalized=False,
                    weight="distance",
                )
                # Average contribution across unordered DEG-gene pairs.
                n_pairs = (len(deg_gene_nodes) * (len(deg_gene_nodes) - 1)) / 2.0
                if n_pairs > 0:
                    for node, val in raw_map.items():
                        metric_map[node] = float(val) / n_pairs
        else:
            metric_map = ga.betweenness_centralities(adjlist, normalized=True)

        tf_metric = {}
        for tf_node in tf_nodes:
            tf = tf_node[1]
            try:
                val = float(metric_map.get(tf_node, 0.0))
            except (TypeError, ValueError):
                val = 0.0
            if val != val:  # NaN guard
                val = 0.0
            tf_metric[tf] = val
        return tf_metric


    def _write_topology_vs_healthy(disorder_name, disorder_adj, healthy_adj, deg_gene_candidates=None):
        outputs = []
        for metric_name, metric_slug in TOPOLOGY_METRIC_SPECS:
            disorder_metrics = _tf_topology_metric_map(
                disorder_adj,
                metric_name,
                deg_gene_candidates=deg_gene_candidates,
            )
            healthy_metrics = _tf_topology_metric_map(
                healthy_adj,
                metric_name,
                deg_gene_candidates=deg_gene_candidates,
            )

            ranked_tfs = sorted(
                disorder_metrics.items(),
                key=lambda kv: (-kv[1], str(kv[0])),
            )[:TOPOLOGY_TOP_N]

            out_txt = (
                f"results/{disorder_name.lower()}_top{TOPOLOGY_TOP_N}_"
                f"tfs_{metric_slug}_vs_healthy.txt"
            )
            with open(out_txt, "w") as f:
                f.write(f"TF\thealthy_{metric_slug}\n")
                for tf, _ in ranked_tfs:
                    f.write(f"{tf}\t{healthy_metrics.get(tf, 0.0):.10g}\n")
            outputs.append({
                "metric": metric_slug,
                "path": out_txt,
                "rows": len(ranked_tfs),
            })
        return outputs


    def _extract_csrf_token(html_text):
        m = re.search(r'name="csrfmiddlewaretoken" value="([^"]+)"', html_text)
        return m.group(1) if m else None


    def _cluereg_query_via_analysis_form(tf_pos, tf_neg, label, signature_mode="tf"):
        errors = []
        tfgene_value = "Gene targeting" if signature_mode == "gene" else "TF targeting"

        for analysis_url in GRAND_ANALYSIS_URLS:
            try:
                # Step 1: fetch form to initialize CSRF cookie + hidden token.
                form_resp = GRAND_SESSION.get(
                    analysis_url,
                    timeout=GRAND_TIMEOUT,
                )
                form_resp.raise_for_status()

                csrf_token = _extract_csrf_token(form_resp.text)
                if not csrf_token:
                    csrf_token = GRAND_SESSION.cookies.get("csrftoken")
                if not csrf_token:
                    raise RuntimeError("missing CSRF token")

                # Step 2: submit CLUEreg form (TF-targeting or Gene-targeting mode).
                form_payload = {
                    "csrfmiddlewaretoken": csrf_token,
                    "contentup": "\n".join(tf_pos),
                    "contentdown": "\n".join(tf_neg),
                    "tfgene": tfgene_value,
                    "ngenes": str(CLUE_TOP_N),
                    "submit": "Submit",
                }

                result_resp = GRAND_SESSION.post(
                    analysis_url,
                    data=form_payload,
                    headers={"Referer": analysis_url},
                    timeout=(GRAND_TIMEOUT[0], 60),
                )
                result_resp.raise_for_status()

                out_html = f"results/cluereg_results_{label}.html"
                with open(out_html, "w") as f:
                    f.write(result_resp.text)
                print(f"saved CLUEreg web result page: {out_html}")
                return True

            except (requests.exceptions.RequestException, RuntimeError) as exc:
                errors.append(f"{analysis_url} :: {exc}")

        if not GRAND_API_STATE["analysis_warning_reported"]:
            print(
                "warning: GRAND analysis-form fallback failed; "
                "unable to retrieve CLUEreg results from website"
            )
            GRAND_API_STATE["analysis_warning_reported"] = True
        print("analysis fallback errors:\n" + "\n".join(errors))
        return False


    def _cluereg_query(tf_pos, tf_neg, label, signature_mode="tf"):
        if not tf_pos or not tf_neg:
            print(
                f"warning: skipping CLUEreg for {label}; "
                "need both non-empty positive and negative lists"
            )
            return

        # Try legacy JSON API first.
        api_ok = False
        if not GRAND_API_STATE["disabled"]:
            payload = {
                "up_regulators": tf_pos,
                "down_regulators": tf_neg,
            }
            try:
                r = _grand_post("/cluereg", payload)
                with open(f"results/cluereg_results_{label}.json", "w") as f:
                    f.write(r.text)
                api_ok = True
            except RuntimeError as exc:
                if not GRAND_API_STATE["cluereg_warning_reported"]:
                    print("warning: legacy GRAND CLUEreg API unavailable; trying /analysis/ form submission fallback")
                    GRAND_API_STATE["cluereg_warning_reported"] = True
                print(str(exc))

        if api_ok:
            return

        # Fallback: submit directly to GRAND analysis web form.
        _cluereg_query_via_analysis_form(tf_pos, tf_neg, label, signature_mode=signature_mode)


    clue_lists = {}
    diff_scores_by_name = {}
    tf_diff_scores_by_name = {}
    gene_diff_scores_by_name = {}

    print(
        "CLUE signature mode:",
        CLUE_SIGNATURE_MODE,
        "| strict signed:",
        CLUE_STRICT_SIGNED,
        "| top_n:",
        CLUE_TOP_N,
    )

    for name in disorders:

        if name not in tf_grn_by_name:
            continue

        print("\nprocessing disorder:", name)

        detf_net = tf_grn_by_name.get(name, {})
        deg_net = deg_grn_by_name.get(name, {})
        deg_gene_candidates = {
            edge[0]
            for edges in deg_net.values()
            for edge in edges
            if edge
        }
        topology_outputs = _write_topology_vs_healthy(
            name,
            detf_net,
            grn,
            deg_gene_candidates=deg_gene_candidates,
        )

        # TF signatures from DETF_GRN
        tf_diff_scores = _differential_targeting_api(
            grn,
            detf_net,
        )
        tf_diff_scores_by_name[name] = tf_diff_scores
        tf_candidates = set(detf_net.keys())
        tf_pos100, tf_neg100 = _top_n_signature(
            tf_diff_scores,
            key_name="TF",
            n=CLUE_TOP_N,
            strict_signed=CLUE_STRICT_SIGNED,
            allowed_values=tf_candidates,
        )

        # Gene signatures from DEG_GRN
        gene_diff_scores = _differential_gene_targeting_local(
            grn,
            deg_net,
        )
        gene_diff_scores_by_name[name] = gene_diff_scores
        gene_candidates = deg_gene_candidates
        gene_pos100, gene_neg100 = _top_n_signature(
            gene_diff_scores,
            key_name="Gene",
            n=CLUE_TOP_N,
            strict_signed=CLUE_STRICT_SIGNED,
            allowed_values=gene_candidates,
        )

        # Keep the mode-selected list for CLUEreg query/output compatibility.
        if CLUE_SIGNATURE_MODE == "gene":
            pos100, neg100 = gene_pos100, gene_neg100
            diff_scores = gene_diff_scores
        else:
            pos100, neg100 = tf_pos100, tf_neg100
            diff_scores = tf_diff_scores

        clue_lists[name] = {
            "positive": pos100,
            "negative": neg100}
        diff_scores_by_name[name] = diff_scores

        _write_cluereg_file(
            pos100,
            f"results/cluereg_{name.lower()}_positive.txt")

        _write_cluereg_file(
            neg100,
            f"results/cluereg_{name.lower()}_negative.txt")

        # Always export both TF and gene ranked signatures per disorder.
        _write_cluereg_file(
            tf_pos100,
            f"results/cluereg_{name.lower()}_tf_positive.txt",
        )
        _write_cluereg_file(
            tf_neg100,
            f"results/cluereg_{name.lower()}_tf_negative.txt",
        )
        _write_cluereg_file(
            gene_pos100,
            f"results/cluereg_{name.lower()}_gene_positive.txt",
        )
        _write_cluereg_file(
            gene_neg100,
            f"results/cluereg_{name.lower()}_gene_negative.txt",
        )

        # Clear-label exports for the most relevant differential-targeting outputs.
        _write_cluereg_file(
            gene_pos100,
            _cluereg_signature_filename(
                name.lower(),
                "genes",
                "most_positive",
                CLUE_TOP_N,
            ),
        )
        _write_cluereg_file(
            gene_neg100,
            _cluereg_signature_filename(
                name.lower(),
                "genes",
                "most_negative",
                CLUE_TOP_N,
            ),
        )

        mode_entity = "genes" if CLUE_SIGNATURE_MODE == "gene" else "tfs"
        _write_cluereg_file(
            pos100,
            _cluereg_signature_filename(
                name.lower(),
                mode_entity,
                "most_positive",
                CLUE_TOP_N,
            ),
        )
        _write_cluereg_file(
            neg100,
            _cluereg_signature_filename(
                name.lower(),
                mode_entity,
                "most_negative",
                CLUE_TOP_N,
            ),
        )

        _cluereg_query(
            pos100,
            neg100,
            name.lower(),
            signature_mode=CLUE_SIGNATURE_MODE,
        )

        print(
            name,
            "CLUEreg mode-selected entries:",
            len(pos100),
            "positive,",
            len(neg100),
            "negative")
        print(
            name,
            "ranked TF signature entries:",
            len(tf_pos100),
            "top,",
            len(tf_neg100),
            "bottom",
        )
        print(
            name,
            "ranked gene signature entries:",
            len(gene_pos100),
            "top,",
            len(gene_neg100),
            "bottom",
        )
        for topo in topology_outputs:
            print(
                name,
                f"topology file ({topo['metric']}):",
                topo["path"],
                f"rows={topo['rows']}",
            )

    _, _, ffl_rows, fbl_rows = gu.detect_regulatory_loops(
        degs,
        tf_grn_by_name,
        graph_label="tf_grn",
    )

    ga.write_csv(
        ffl_rows,
        "results/deggrn_feedforward_loops.csv",
    )
    ga.write_csv(
        fbl_rows,
        "results/deggrn_feedback_loops.csv",
    )


def test_merge_adjlist(brainother, brainbg, braincerebellum, brains):
    def _edge_set(adj):
        return {(tf, gene) for tf, edges in adj.items() for gene, _ in edges}

    def _edge_weights(adj):
        return {(tf, gene): w for tf, edges in adj.items() for gene, w in edges}

    edges_other = _edge_set(brainother)
    edges_bg = _edge_set(brainbg)
    edges_cerebellum = _edge_set(braincerebellum)
    edges_merged = _edge_set(brains)

    # 1) merged contains all edges from all inputs
    assert edges_other.issubset(edges_merged)
    assert edges_bg.issubset(edges_merged)
    assert edges_cerebellum.issubset(edges_merged)
    assert edges_merged == edges_other | edges_bg | edges_cerebellum

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
