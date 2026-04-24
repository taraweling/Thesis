import matplotlib.pyplot as plt
from pyvis.network import Network
import math
import graph_utils as gu

# Shared process-level cache for ENSG -> display label lookups.
_GENE_LABEL_CACHE = {}

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

def pyviz_deggrn(adjlist, outfile="grn.html", directed=True, top_tfs=None, top_degs=None): # option to make not directed so it resembles a co expression net
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
          TF nodes     > triangle
          Gene nodes   > circle
          Node size    >  log10(pval)
          Node color   > log2fc (red up / blue down)
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

            node_info.setdefault(gene, {
                "type": "Gene",
                "log2fc": log2fc,
                "pval": pval
            })

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
        color = "#cccccc"
        if info.get("log2fc") is not None:
            if info["log2fc"] > 0:
                color = "#d73027"
            elif info["log2fc"] < 0:
                color = "#4575b4"

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
            weight = abs(edge[1]) if len(edge) > 1 else 1.0
            G.add_edge(tf, gene, value=weight)

    # Write output
    if len(G.edges) == 0:
        print("Warning: pyviz_deggrn has no edges to render.")
        return

    G.toggle_physics(True)
    G.show_buttons(filter_=["physics"])
    G.write_html(outfile)
    return 

def viz_graph(edgelist, outfile, top_tfs=None, top_degs=None): # constructed using anna ritz's course assignment as inspiration
    """
    Visualize a directed graph from an enriched edgelist and write it to an HTML file.

    Input rows:
        [TF, Gene, weight, disorder, study, year, tissue, log2fc, pval]

    Visualization (if present):
          Edge width   > |weight|
          TF nodes     > triangle
          Gene nodes   > circle
          Node size    >  log10(pval)
          Node color   > log2fc (red up / blue down)
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
        node_info.setdefault(gene, {
            "type": "Gene",
            "log2fc": log2fc,
            "pval": pval})

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
        color = "#cccccc"
        if info.get("log2fc") is not None:
            if info["log2fc"] > 0:
                color = "#d73027"   # upregulated
            elif info["log2fc"] < 0:
                color = "#4575b4"   # downregulated

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
        weight = abs(row[2]) if len(row) > 2 else 1.0

        G.add_edge(tf, gene, value=weight)

    # Output
    if len(G.edges) == 0:
        print("Warning: viz_graph has no edges to render.")
        return
    G.toggle_physics(True)
    G.show_buttons(filter_=["physics"])
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
