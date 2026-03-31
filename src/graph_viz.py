import matplotlib.pyplot as plt
from pyvis.network import Network
import math

# All functions now assume edgelist rows look like:
# [TF, Gene, weight] OR
# [TF, Gene, weight, disorder, study, year, tissue, log2fc, pval]

# Visualize
## MUST FIGURE OUT A WAY TO VISUALIZE ONLY THE GENES THAT THE MOST TFS POINT TO 

#  returns graph of the most differentally expressed TFs 
def regulon_viz(adjlist, outfile ="regulon.html", top_tfs=20, top_genes_by_tfs=None, directed=True): 
    # rank tfs by mean |log2fc| of their targets
    adjlist = _top_tfs_by_abs_log2fc_adjlist(adjlist, top_tfs)
    # keep genes hit by the most tfs
    adjlist = _top_genes_by_tfs_adjlist(adjlist, top_genes_by_tfs)
    return pyviz_deggrn(adjlist, outfile=outfile, directed=directed)
    
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

def _top_tfs_by_abs_log2fc_adjlist(adjlist, top_tfs):
    if top_tfs is None:
        return adjlist

    tf_scores = {}
    for tf, edges in adjlist.items():
        fcs = []
        targets = set()
        for edge in edges:
            if isinstance(edge, (list, tuple)):
                targets.add(edge[0])
                if len(edge) > 7 and edge[7] is not None:
                    try:
                        fcs.append(abs(float(edge[7])))
                    except (TypeError, ValueError):
                        pass
            else:
                targets.add(edge)
        if fcs:
            tf_scores[tf] = sum(fcs) / len(fcs)
        else:
            tf_scores[tf] = len(targets)

    top = sorted(tf_scores, key=tf_scores.get, reverse=True)[:top_tfs]
    return {tf: adjlist[tf] for tf in top}

def _top_genes_by_tfs_adjlist(adjlist, top_genes):
    if top_genes is None:
        return adjlist

    gene_to_tfs = {}
    for tf, edges in adjlist.items():
        for edge in edges:
            gene = edge[0] if isinstance(edge, (list, tuple)) else edge
            gene_to_tfs.setdefault(gene, set()).add(tf)

    top = sorted(gene_to_tfs, key=lambda g: len(gene_to_tfs[g]), reverse=True)[:top_genes]
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

def _top_degs_adjlist(adjlist, top_degs):
    if top_degs is None:
        return adjlist

    # Rank genes by |log2fc| desc, then pval asc (best evidence across edges)
    gene_best = {}
    for tf, edges in adjlist.items():
        for edge in edges:
            gene = edge[0] if isinstance(edge, (list, tuple)) else edge
            log2fc = edge[7] if isinstance(edge, (list, tuple)) and len(edge) > 7 else None
            pval = edge[8] if isinstance(edge, (list, tuple)) and len(edge) > 8 else None

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

def _top_genes_by_tfs_edgelist(edgelist, top_genes):
    if top_genes is None:
        return edgelist

    gene_to_tfs = {}
    for row in edgelist:
        tf = row[0]
        gene = row[1]
        gene_to_tfs.setdefault(gene, set()).add(tf)

    top = sorted(gene_to_tfs, key=lambda g: len(gene_to_tfs[g]), reverse=True)[:top_genes]
    keep = set(top)
    return [row for row in edgelist if row[1] in keep]

def pyviz_deggrn(adjlist, outfile="grn.html", directed=True, top_tfs=None, top_degs=None, top_genes_by_tfs=None): # option to make not directed so it resembles a co expression net
    """
    Visualize a Gene Regulatory Network (GRN) stored as an adjacency list.

    Input formats supported:

    1) Simple:
        {
            TF: [(Gene, weight), ...]
        }

    2) Annotated:
        {
            TF: [
                (Gene, weight, disorder, study, year, tissue, log2fc, pval),
                ...
            ]
        }

    Output:
        Interactive HTML visualization written to `outfile`.

    Visual encodings (if available):
          Edge width   > |weight|
          TF nodes     > triangle
          Gene nodes   > circle
          Node size    >  log10(pval)
          Node color   > log2fc (red up / blue down)
    Filtering:
          top_tfs  > keep TFs with most unique targets
          top_degs > keep genes with largest |log2fc|, then lowest pval
          top_genes_by_tfs > keep genes with most TF regulators
    """

    adjlist = _top_tfs_adjlist(adjlist, top_tfs)
    adjlist = _top_degs_adjlist(adjlist, top_degs)
    adjlist = _top_genes_by_tfs_adjlist(adjlist, top_genes_by_tfs)
    G = Network(directed=directed)

    # Collect node metadata

    node_info = {}

    for tf, targets in adjlist.items():
        node_info.setdefault(tf, {"type": "TF"})

        for edge in targets:
            gene = edge[0]
            weight = edge[1] if len(edge) > 1 else None

            log2fc = edge[7] if len(edge) > 7 else None
            pval = edge[8] if len(edge) > 8 else None

            node_info.setdefault(gene, {
                "type": "Gene",
                "log2fc": log2fc,
                "pval": pval
            })

    # Add nodes    

    for node, info in node_info.items():

        shape = "triangle" if info.get("type") == "TF" else "circle"

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
            label=node,
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

def viz_graph(edgelist, outfile, top_tfs=None, top_degs=None, top_genes_by_tfs=None): # constructed using anna ritz's course assignment as inspiration
    """
    Visualize a directed graph from an enriched edgelist and write it to an HTML file.

    Input rows:
        [TF, Gene, weight, disorder, study, year, tissue, log2fc, pval]

    Visualization (if present):
          Edge width   > |weight|
          TF nodes     > triangle
          Gene nodes   > circle
          Node size    >  log10(pval)
          Node outline color  > disorder
          Node fill color  > log2fc (red up / blue down)
    Filtering:
          top_tfs  > keep TFs with most targets
          top_degs > keep genes with largest |log2fc|, then lowest pval
          top_genes_by_tfs > keep genes with most TF regulators
    """

    if top_tfs is not None:
        tf_counts = {}
        for row in edgelist:
            tf = row[0]
            tf_counts[tf] = tf_counts.get(tf, 0) + 1
        top = set(sorted(tf_counts, key=tf_counts.get, reverse=True)[:top_tfs])
        edgelist = [row for row in edgelist if row[0] in top]
    edgelist = _top_degs_edgelist(edgelist, top_degs)
    edgelist = _top_genes_by_tfs_edgelist(edgelist, top_genes_by_tfs)

    G = Network(directed=True)

    # Collect nodes + metadata
    node_info = {}

    for row in edgelist:
        tf = row[0]
        gene = row[1]

        weight = row[2] if len(row) > 2 else None
        disorder = row[3] if len(row) > 3 else None
        log2fc = row[7] if len(row) > 7 else None
        pval = row[8] if len(row) > 8 else None

        # TF node
        node_info.setdefault(tf, {
            "type": "TF",
            "disorder": disorder,
            "log2fc": log2fc,
            "pval": pval})

        # Gene node
        node_info.setdefault(gene, {
            "type": "Gene",
            "disorder": disorder,
            "log2fc": log2fc,
            "pval": pval})

    # Disorder → outline color map
    disorders = sorted({
        info["disorder"]
        for info in node_info.values()
        if info["disorder"] is not None})

    palette = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33","#a65628", "#f781bf"]

    disorder_color = {d: palette[i % len(palette)] for i, d in enumerate(disorders)}

    # Add nodes with styling

    for n, info in node_info.items():

        # shape
        shape = "triangle" if info["type"] == "TF" else "circle"

        # size from pval
        size = 20
        if info["pval"] is not None and info["pval"] > 0:
            import math
            size = 10 + (-math.log10(info["pval"]) * 6)
            size = max(8, min(60, size))

        # fill color from log2fc
        fill = "#cccccc"
        if info["log2fc"] is not None:
            if info["log2fc"] > 0:
                fill = "#d73027"   # upregulated
            elif info["log2fc"] < 0:
                fill = "#4575b4"   # downregulated

        # outline color from disorder
        border = disorder_color.get(info["disorder"], "#000000")

        G.add_node(
            n,
            label=n,
            shape=shape,
            size=size,
            color={
                "background": fill,
                "border": border
            },
            borderWidth=3
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
