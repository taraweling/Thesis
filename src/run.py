import graph_utils as gu
import graph_viz as gv
import graph_algos as ga

RUN_TESTS = False
TOP_DEGS = 10
# If True, override TOP_DEGS per network with its unique gene count.
TOP_DEGS_AUTO = True

# INSTRUCTIONS: RUN WHILE IN SRC FOLDER
def _unique_gene_count(adjlist):
    genes = set()
    for edges in adjlist.values():
        for edge in edges:
            gene = edge[0] if isinstance(edge, (list, tuple)) else edge
            genes.add(gene)
    return len(genes)

def _resolve_top_degs(adjlist):
    if TOP_DEGS_AUTO:
        return _unique_gene_count(adjlist)
    return TOP_DEGS

def main():

    # obtain adjlist for each disorder (handle in main or helper fn?)
    ## what information do I want in my graph?
    """
    as a bipartite graph, node attribute should either be TF or gene.
    the TFs are DEGs and the target genes can be either DEGs or not differentially expressed genes in the GRAND regulatory network
    """

    # get the combined GRAND GRN (what check can I use?)
    brainother = gu.make_adjlist('data/Brain_Other.csv', 1.7)
    brainbg = gu.make_adjlist('data/Brain_Basal_Ganglia.csv', 1.7)
    braincb = gu.make_adjlist('data/Brain_Cerebellum.csv',1.7)
    
    # merge the two brain GRNs and test the merge behavior
    brain = gu.merge_adjlist(brainother, brainbg)
    brains = gu.merge_adjlist(brain,braincb)
    
    if RUN_TESTS:
        test_merge_adjlist(brainother, brainbg, brain)
        test_merge_adjlist(brain, braincb, brains)
        
    grn = gu.ensemblify(brains)
    
    # run.py (right after ensemblify)
    ensg_tfs = [tf for tf in grn if gu.ENSEMBL_ID_RE.match(tf)]
    print("TFs mapped to ENSG:", len(ensg_tfs), "of", len(grn))

    # get adjlists per disorder by inputting the name of the disorder
    ## (options = AD, ADHD, BD, SZ, MDD, OCD) and the file location, outputting a list of lists
    disorders = ['AD', 'ADHD', 'ASD', 'BD', 'MDD', 'OCD', 'SZ']
    data_path = 'data/DEGDataStrictLFC.csv'
    # can be changed to DEGData
    degs = {d: gu.disorder_list(data_path, d) for d in disorders}

    # per-disorder outputs (the by-name suffix allows us to 
    detf_deggrn_by_name = {}
    tf_grn_by_name = {}
    deg_grn_by_name = {}
    detf_deggrn_edgelist_by_name = {}
    tf_grn_edgelist_by_name = {}
    reg_scores_by_name = {}
    tf_regulators_by_name = {}
    edgeweight_summary_by_name = {}
    log2fc_summary_by_name = {}
    summary_rows = []
    overlap_rows = []
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
        overlap_rows.append(
            gu.overlap_summary_row(name, degset, grn_tfs, grn_nodes))

        # build DEG-filtered GRNs per disorder
        tf_grn = gu.de_grn_tfsonly(grn, data)    # DETFs only
        deg_grn = gu.de_grn_degsonly(grn, data)  # DEG targets only

        # Only build detf_deggrn if it would differ from tf_grn.
        detf_is_tf = True
        for tf, edges in tf_grn.items():
            for edge in edges:
                gene = edge[0] if isinstance(edge, (list, tuple)) else edge
                if gene not in degset:
                    detf_is_tf = False
                    break
            if not detf_is_tf:
                break

        if detf_is_tf:
            detf_deggrn = tf_grn
        else:
            detf_deggrn = gu.de_grn_both(grn, data)  # DETFs and DEG targets

        # store by disorder name
        detf_deggrn_by_name[name] = detf_deggrn
        tf_grn_by_name[name] = tf_grn
        deg_grn_by_name[name] = deg_grn

        # naming system for quick access (e.g., adhd_detf_deggrn)
        globals()[f"{name.lower()}_detf_deggrn"] = detf_deggrn
        globals()[f"{name.lower()}_tf_grn"] = tf_grn
        globals()[f"{name.lower()}_deg_grn"] = deg_grn

        # apply graph_algos to each disorder-specific GRN
        reg_scores_by_name[name] = ga.regulatory_scores(detf_deggrn)
        detf_deggrn_edgelist_by_name[name] = gu.adjlist2edgelist(detf_deggrn)
        tf_grn_edgelist_by_name[name] = gu.adjlist2edgelist(tf_grn)
        tf_regulators_by_name[name] = ga.regulator_detection(tf_grn, data)
        edgeweight_summary_by_name[name] = ga.edgeweight_summary(grn, data)
        log2fc_summary_by_name[name] = {
            "detf_deggrn": ga.log2fc_summary(detf_deggrn),
            "tf_grn": ga.log2fc_summary(tf_grn),
            "deg_grn": ga.log2fc_summary(deg_grn),}

        # collect per-disorder summary row for CSV
        summary_rows.append(
            ga.disorder_summary_row(
                name=name,
                disorderlist=data,
                detf_deggrn=detf_deggrn,
                tf_grn=tf_grn,
                deg_grn=deg_grn,
                reg_scores=reg_scores_by_name[name],
                tf_regulators=tf_regulators_by_name[name],
                edgeweight_summary=edgeweight_summary_by_name[name],))

        # sets for pairwise Jaccard comparisons
        deg_gene_sets[name] = {rec[0] for rec in data}
        detf_edge_sets[name] = {
            (tf, edge[0])
            for tf, edges in detf_deggrn.items()
            for edge in edges}
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
                for r in tf_scores[:n])

        print(f"\n{name} summary stats:")
        if detf_is_tf:
            print("  detf_deggrn is identical to tf_grn; skipping separate build.")
        else:
            print(f"  detf_deggrn TFs={len(detf_deggrn)} edges={_edge_count(detf_deggrn)}")
        print(f"  tf_grn TFs={len(tf_grn)} edges={_edge_count(tf_grn)}")
        print(f"  deg_grn TFs={len(deg_grn)} edges={_edge_count(deg_grn)}")
        print(f"  top DETF regulators: {_top_regs(reg_scores_by_name[name])}")
        print(f"  top TF regulators: {_top_regs(tf_regulators_by_name[name])}")
        print(f"  edgeweight summary: {edgeweight_summary_by_name[name]}")
        print(f"  log2fc summary (detf_deggrn): {log2fc_summary_by_name[name]['detf_deggrn']}")
        print(f"  log2fc summary (tf_grn): {log2fc_summary_by_name[name]['tf_grn']}")
        print(f"  log2fc summary (deg_grn): {log2fc_summary_by_name[name]['deg_grn']}")

    # visualize graphs per disorder
        if not detf_is_tf:
            detf_deggrnedgelist = detf_deggrn_edgelist_by_name[name]
            top_degs_detf = _resolve_top_degs(detf_deggrn)
            gv.viz_graph(
                detf_deggrnedgelist,
                f"results/deggrn_{name.lower()}.html",
                top_degs=top_degs_detf,)

        # addition of deggrn pyviz rather than detf_deggrn
        deg_grnedgelist = deg_grn_by_name[name]
        top_degs_deg = _resolve_top_degs(deg_grn)
        gv.pyviz_deggrn(
            deg_grn,
            outfile=f"results/pyviz_{name.lower()}_deggrn.html",
            top_degs=top_degs_deg,)
        
        if not detf_is_tf:
            gv.pyviz_deggrn(
                detf_deggrn,
                outfile=f"results/pyviz_{name.lower()}.html",
                top_degs=top_degs_detf,)

    # write CSV comparisons
    ga.write_csv(summary_rows, "results/deggrn_disorder_summary.csv")
    gu.write_overlap_summary(
        overlap_rows,
        "results/deggrn_overlap_summary.csv",)

    # overlay network across disorders (TF-only GRNs)
    overlay_tf_edges = []
    for edgelist in tf_grn_edgelist_by_name.values():
        overlay_tf_edges.extend(edgelist)
    gv.viz_overlay_disorders(
        overlay_tf_edges,
        "results/overlay_tf_grn_disorders.html",
        top_degs=None,) # changed from top_degs=TOP_DEGS
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
        assert abs(wm - ((w1 + w2) / 2)) < 1e-9 # check thqt 


if __name__ == '__main__':
    main()
