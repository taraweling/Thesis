import matplotlib.pyplot as plt
from pyvis.network import Network
import networkx as nx
import numpy as np
import graph_algos as ga
import graph_utils as gu
import json
import py4cytoscape as pyc

# All functions now assume edgelist rows look like:
# [TF, Gene, weight] OR
# [TF, Gene, weight, disorder, study, year, tissue, log2fc, pval]


# Visualize
## MUST FIGURE OUT A WAY TO VISUALIZE ONLY THE GENES THAT THE MOST TFS POINT TO 
from pyvis.network import Network

def lastditch(adjlist, outfile="grn.html", directed=True):
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
        - Edge width  -> |weight|
        - TF nodes    -> triangle
        - Gene nodes  -> circle
        - Node size   -> -log10(pval)
        - Node color  -> log2fc (red up / blue down)
    """

    G = Network(directed=directed)

    # ------------------------------------------------------------------
    # Collect node metadata
    # ------------------------------------------------------------------

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

    # ------------------------------------------------------------------
    # Add nodes
    # ------------------------------------------------------------------

    import math

    for node, info in node_info.items():

        shape = "triangle" if info.get("type") == "TF" else "circle"

        # size from p-value
        size = 20
        if info.get("pval") is not None and info["pval"] > 0:
            size = min(60, 10 + (-math.log10(info["pval"]) * 6))

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

    # ------------------------------------------------------------------
    # Add edges
    # ------------------------------------------------------------------

    for tf, targets in adjlist.items():
        for edge in targets:
            gene = edge[0]
            weight = abs(edge[1]) if len(edge) > 1 else 1.0
            G.add_edge(tf, gene, value=weight)

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------

    G.toggle_physics(True)
    G.show_buttons(filter_=["physics"])
    G.write_html(outfile)


def visualize_deg_grn(adjlist, min_weight=None, max_pval=None, min_abs_log2fc=None, max_edges=1000, seed=0):
    """
    Visualize a DEG-GRN adjacency list.

    adjlist format:
    {
        TF: [
            (Gene, GRN_weight, disorder, study, year, tissue, log2fc, pval),
            ...
        ]
    }

    Filters:
        min_weight      → drop weak GRN edges
        max_pval        → drop nonsignificant DEGs
        min_abs_log2fc  → drop weak expression changes
        max_edges       → subsample edges if still huge
    """

    G = nx.DiGraph()

    # Build graph with filters

    for tf, edges in adjlist.items():
        for gene, w, dis, study, year, tissue, lfc, pval in edges:

            if min_weight is not None and abs(w) < min_weight:
                continue

            if max_pval is not None and pval > max_pval:
                continue

            if min_abs_log2fc is not None and abs(lfc) < min_abs_log2fc:
                continue

            G.add_edge(
                tf,
                gene,
                weight=w,
                disorder=dis,
                study=study,
                year=year,
                tissue=tissue,
                log2fc=lfc,
                pval=pval,
            )

    # Downsample if enormous

    if G.number_of_edges() > max_edges:
        import random
        random.seed(seed)

        edges = list(G.edges(data=True))
        keep = random.sample(edges, max_edges)

        G2 = nx.DiGraph()
        for u, v, d in keep:
            G2.add_edge(u, v, **d)

        G = G2

    # Layout

    pos = nx.spring_layout(G, seed=seed, k=1.2)

    # TF nodes = sources
    tfs = set(adjlist.keys())
    tf_nodes = [n for n in G.nodes if n in tfs]
    gene_nodes = [n for n in G.nodes if n not in tfs]

    # Draw

    plt.figure(figsize=(14, 12))

    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=tf_nodes,
        node_color="firebrick",
        node_size=800,
        label="TFs",
    )

    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=gene_nodes,
        node_color="steelblue",
        node_size=400,
        label="Targets",
    )

    nx.draw_networkx_edges(
        G,
        pos,
        arrows=True,
        alpha=0.5,
        width=1,
    )

    nx.draw_networkx_labels(G, pos, font_size=7)

    plt.title(
        f"DEG–GRN network\nNodes={G.number_of_nodes()}  Edges={G.number_of_edges()}"
    )

    #plt.legend()
    plt.axis("off")
    plt.tight_layout()
    plt.savefig('results/deggrn.jpeg')
    plt.show()

def viz_graph(edgelist, outfile): # constructed using anna ritz's course assignment as inspiration
    """
    Visualize a directed graph from an enriched edgelist and write it to an HTML file.

    Input rows:
        [TF, Gene, weight, disorder, study, year, tissue, log2fc, pval]

    Visualization (if present):
        - Edge width  -> |weight|
        - TF nodes    -> triangle
        - Gene nodes  -> circle
        - Node size   -> -log10(pval)
        - Node outline color -> disorder
        - Node fill color -> log2fc (red up / blue down)
    """

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
            size = min(60, 10 + (-math.log10(info["pval"]) * 6))

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
    G.toggle_physics(True)
    G.show_buttons(filter_=["physics"])
    G.write_html(outfile)

def old_viz_graph(edgelist, outfile):
    """
    Visualize a directed graph represented as an enriched edgelist
    and write it to an HTML file.

    Input rows:
        [TF, Gene, weight, ...metadata...]

    Only TF -> Gene is used for topology.
    Weight is optionally used to scale edge width.
    """

    G = Network(directed=True)

    nodes = set()
    for row in edgelist:
        u = row[0]
        v = row[1]
        nodes.add(u)
        nodes.add(v)

    for n in nodes:
        G.add_node(n, label=n, color="#FFFFFF")

    for row in edgelist:
        u = row[0]
        v = row[1]
        weight = row[2] if len(row) >= 3 else 1.0
        G.add_edge(u, v, value=abs(weight))

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
