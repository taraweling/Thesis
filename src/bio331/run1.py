from graph_class import process_data
import graph_vis as gv 
import graph_algos as ga 
import p4utils as utils
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex # colors for graph
import numpy as np
from scipy.sparse.linalg import eigsh
from scipy.sparse import csr_matrix
from scipy.sparse import csgraph
import scipy.stats
from scipy.stats import mannwhitneyu
from pyvis.network import Network
import sys

# by tara

def main():
    
    # code to obtain adjlist, inputs location of csvs and threshold of distances
    e_adjlist = process_data('adjlistdata/epithelial.csv',1) # 10 most epithelial slices
    m_adjlist = process_data('adjlistdata/mesenchymal.csv',1) # 10 most mesenchymal slices
    
    # iris's code graph_vis calls thomas's code ga
    
    outfile1 = 'boxplots.png'
    outfile2 = 'violinplots.png'
    
    # code to visualize graph, calls graph_viz.py and writes to outfile 1, 2
    deg_data, clust_data, close_data, eccen_data, betw_data = gv.net_compare(e_adjlist, m_adjlist, outfile1, outfile2) 
    
    #print('\n\n',m_eccentricity,'\n\n',e_eccentricity)
    outfile = 'adjlistdata/'
    
    # turn below code into a function using dictionaries to pair file name with title
    epi1, mes1 = measure2list(betw_data)
    compare_dist(epi1,mes1,'Betweenness',outfile)
    
    epi2, mes2 = measure2list(deg_data)
    compare_dist(epi2,mes2,'Degree',outfile)
    
    epi3, mes3 = measure2list(clust_data)
    compare_dist(epi3,mes3,'Clustering_Coefficient',outfile)
    
    epi4, mes4 = measure2list(close_data)
    compare_dist(epi4,mes4,'Closeness_Centrality',outfile)
    
    epi5, mes5 = measure2list(eccen_data)
    compare_dist(epi5,mes5,'Eccentricity_Centrality',outfile)
    
    #### visualize graph (girvan newman clustering seemed to fail)
    viz_graph(e_adjlist,'epithelial_graph.html')
    viz_graph(m_adjlist,'mesenchymal_graph.html')
    
    #### higher order clustering doesn't work
    # e_lap = laplacian_matrix(e_adjlist,"Epithelial")
    # m_lap = laplacian_matrix(m_adjlist,"Mesenchymal")
    
    #### conduct statistical testing 
    # cliff's delta = magnitude and direction of EMT associated change (>=0.47 = large difference)
    
    print("betweenness")
    results1 = stat_test(epi1,mes1)
    print("U statistic:", results1['U_statistic'])
    print("p-value:", results1['p_value'])
    print("Cliff's delta:", results1['cliffs_delta'])
    print("\n")
    
    print("degree distribution")
    results2 = stat_test(epi2,mes2)
    print("U statistic:", results2['U_statistic'])
    print("p-value:", results2['p_value'])
    print("Cliff's delta:", results2['cliffs_delta'])
    print("\n")
    
    print("clustering coefficient")
    results3 = stat_test(epi3,mes3)
    print("U statistic:", results3['U_statistic'])
    print("p-value:", results3['p_value'])
    print("Cliff's delta:", results3['cliffs_delta'])
    print("\n")
    
    print("closeness centrality")
    results4 = stat_test(epi4,mes4)
    print("U statistic:", results4['U_statistic'])
    print("p-value:", results4['p_value'])
    print("Cliff's delta:", results4['cliffs_delta'])
    print("\n")
    
    print("eccentricity centrality")
    results4 = stat_test(epi4,mes4)
    print("U statistic:", results4['U_statistic'])
    print("p-value:", results4['p_value'])
    print("Cliff's delta:", results4['cliffs_delta'])
    print("\n")
    return

def measure2list(measure): # outputs of gv.net_compare better be ordered
    epi = [measure[0]]
    mes = [measure[1]]
    [epi] = epi
    [mes] = mes
    epiavg = sum(epi) / len(epi)
    mesavg = sum(mes) / len(mes)
    return epi, mes

def dist(measure):
    # Calculate the frequency of each measure
    degree_counts = np.bincount(measure)
    total_nodes = len(measure)
    
    # Normalize the counts to get probabilities, necessary bc epithelial and mesenchymal adjlists are different sizes
    degree_probabilities = degree_counts / total_nodes
    
    return degree_probabilities

def compare_dist(epi, mes, value, outfile):
    # Calculate normalized measure distributions
    prob1 = dist(epi)
    prob2 = dist(mes)
    
    # Plot the measure's distributions
    plt.figure(figsize=(8, 6))
    plt.plot(prob1, label='Epithelial', color='red')
    plt.plot(prob2, label='Mesenchymal', color='blue')
    
    # turn below into its own fn
    plt.xlabel(value)
    plt.ylabel('Probability')
    plt.legend()
    tt = value + ' Distribution Comparison'
    plt.title(tt)
    
    save = outfile + value
    plt.savefig(str(save))
    return 

def get_clusters(adjlist): # girvan-newman clustering from p4
    nodes = ga.get_nodes(adjlist)
    edgelist = []
    for node, neighbors in adjlist.items():
        for neighbor in neighbors:
            edgelist.append((node, neighbor)) # ensure edgelist is tuple
    partitions = [[nodes]]  # start with all nodes in a single cluster
    
    while len(edgelist) > 0:
        # compute edge betweenness
        paths = utils.all_pairs_shortest_paths(adjlist)
        edge_bw = utils.edge_betweenness(nodes, edgelist, paths)
        
        # C3: remove edge with highest betweenness
        highest = max(edge_bw, key=edge_bw.get)
        forward = list(highest)
        reverse = list(highest)
        reverse.reverse()
        
        utils.remove_from_adjlist(forward, adjlist)
        utils.remove_from_edgelist(forward, edgelist)
        utils.remove_from_adjlist(reverse, adjlist)
        
        # C4: check if removal split a component
        list1 = utils.conn_comp(adjlist, forward[0])
        list2 = utils.conn_comp(adjlist, forward[1])
        
        if list1 != list2:
            old_cluster = partitions.pop()  # remove the cluster we just split
            partitions.extend(utils.split_partition(old_cluster, list1, list2))

    return partitions

def viz_graph(weighted_adjlist,outfile): # inputs adjlist, outputs graph
    # build graph
    G = Network()
    """for node, neighbors in weighted_adjlist.items():
        if node not in G.get_nodes():
            G.add_node(node)
        for neighbor, weight in neighbors:
            if neighbor not in G.get_nodes():
                G.add_node(neighbor)
            G.add_edge(node, neighbor, weight=weight)"""
    
    all_nodes = set(weighted_adjlist.keys()) # add all nodes
    for neighbors in weighted_adjlist.values():
        for neighbor, _ in neighbors:
            all_nodes.add(neighbor)

    for node in all_nodes:
        G.add_node(str(node))
    
    # add all edges
    for node, neighbors in weighted_adjlist.items():
        for neighbor, weight in neighbors:
            G.add_edge(str(node), str(neighbor), weight=weight)

    nodes = ga.get_nodes(weighted_adjlist)
    
    adjlist_unweighted = ga.unweight_adjlist(weighted_adjlist,nodes)
    # draw graph
    G.save_graph(outfile)
    """ 
    # detect clusters
    clusters = get_clusters(adjlist_unweighted)
    
    # Assign cluster colors
    colormap = plt.cm.get_cmap('tab10', len(clusters))
    cluster_colors = {}
    for i, cluster in enumerate(clusters):
        for node in cluster:
            cluster_colors[str(node)] = to_hex(colormap(i))
    
    # Update nodes in PyVis
    for node in G.nodes:
        node_id = str(node['id'])
        if node_id in cluster_colors:
            node['color'] = cluster_colors[node_id]
        else:
            node['color'] = '#cccccc'  # default color for missing nodes

        # Optional: size nodes by degree
        node_degree = 0
        for e in G.edges:
            if e['from'] == node_id or e['to'] == node_id:
                node_degree += 1

        node['size'] = 10 + node_degree * 2  # base size + scaled by degree so important nodes stand out
    """
    
    return G

def stat_test(epi, mes, ): # compare node degree dists bw epi and mes using mann-whitney u test
    # uses numpy = convert to numpy arrays
    epi = np.asarray(epi, dtype=float)
    mes = np.asarray(mes, dtype=float)
    
    # removes NaNs
    epi = epi[~np.isnan(epi)]
    mes = mes[~np.isnan(mes)]
    
    # runs Mann-Whitney U test (non-parametric)
    U, p_val = mannwhitneyu(epi,mes,alternative='two-sided') # replace with 'less' or 'greater' for alt hypothesis
    
    greater = 0
    less = 0
    
    # compute effect size = cliff's delta
    for e in epi:
        greater += np.sum(e > mes)
        less += np.sum(e < mes)

    cliff_delta = (greater - less) / (len(epi) * len(mes))

    results = {
        'U_statistic': U,
        'p_value': p_val,
        'cliffs_delta': cliff_delta,
        'epi_median_degree': np.median(epi),
        'mes_median_degree': np.median(mes),
        'epi_mean_degree': np.mean(epi),
        'mes_mean_degree': np.mean(mes),
        'n_epi': len(epi),
        'n_mes': len(mes)
    }
    # returns dictionary of test stat, p-value, effect size and summary stats
    return results
# unused currently taken from https://matplotlib.org/stable/gallery/lines_bars_and_markers/scatter_hist.html + edited then debugged with chatgpt bc idk matplotlib
def scatter_hist(x, y, ax, ax_histx, ax_histy):
    # No labels on histograms
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # Scatter plot
    ax.scatter(x, y)

    # Determine nice limits for histograms by hand
    binwidth = 0.25
    xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    lim = (int(xymax / binwidth) + 1) * binwidth

    # Create bins for the histograms
    bins = np.arange(-lim, lim + binwidth, binwidth)
    
    # Plot the histograms on the respective axes
    ax_histx.hist(x, bins=bins)
    ax_histy.hist(y, bins=bins, orientation='horizontal')
    
    return 

# below doesn't work bc laplacian is full of infs
def laplacian_matrix(adjlist,state):
    # first convert adjmat to lapmat
    adjmat = ga.adjlist_to_adjmat(adjlist)
    adjmat = np.asarray(adjmat, dtype=float)
    
    degrees = np.sum(adjmat, axis=1)
    laplacian_matrix = np.diag(degrees) - adjmat
    # degree_matrix = np.diag(degrees)
    # laplacian_matrix = degree_matrix - adjmat

    # get eigenvalues and eigenvectors by doing something called eigendecomposition
    eigenvalues, eigenvectors = np.linalg.eigh(laplacian_matrix)
    # 'eigenvectors' is the eigenmatrix (column v[:, i] is the eigenvector corresponding to the eigenvalue w[i])
    # 'eigenvalues' are the corresponding eigenvalues

    # sanity check - thank you chatgpt
    n = laplacian_matrix.shape[0]
    print("Laplacian size:", n)
    #print("Laplacian shape:", laplacian_matrix.shape)
    print(type(laplacian_matrix))
    # --- remove isolated nodes (degree == 0) ---
    keep = degrees > 0
    if keep.sum() < 4:
        raise ValueError("Graph too small after removing isolated nodes")
    L = laplacian_matrix[keep][:, keep]
    # --- convert to sparse ---
    L = csr_matrix(L)
    # --- enforce symmetry (important!) ---
    L = 0.5 * (L + L.T)
    # find fiedler's value = second smallest eigenvalue --> functionally the  minimum
    L = largest_connected_component(L)
    print(L)
    #eigenvalues2 = eigsh(laplacian_matrix, k=2, which='SM', return_eigenvectors=False)
    # eigenvalues2 = eigsh(laplacian_matrix, k=2, which='SA')
    # eigenvalues2 = np.sort(eigenvalues2)
    # fiedler_value = eigenvalues2[1] # second lowest eigenvalue
    # lower_bound = fiedler_value / 2
    # upper_bound = np.sqrt(2 * fiedler_value)
    # --- compute Fiedler value safely ---
    #eigenvalues2 = eigsh(L,k=2,sigma=1e-8,which="LM",return_eigenvectors=False) # shift slightly away from zero.
    eigenvalues2 = np.linalg.eigvalsh(L)

    eigenvalues2 = np.sort(eigenvalues2)
    fiedler_value = eigenvalues2[1]
    lower_bound = fiedler_value / 2
    upper_bound = np.sqrt(2 * fiedler_value)
    print("Fiedler's value: ",fiedler_value,lower_bound,upper_bound)

    ## find spectral radius by finding the maximum |eigenvalue|. pseudocode below
    #sr = max(abs(eigenvalues, key=abs))
    # --- spectral radius (dense, only if needed) ---
    eigenvalues_full = np.linalg.eigvalsh(L.toarray())
    sr = np.max(np.abs(eigenvalues_full))
    
    ### if sr < 1, system is stable (technically something called schur stable 
    if sr < 1:
        print(state," System is Schur stable: ", sr)
    else:
        print(state, " System is not Schur stable: ",sr)

    # from here we can print out the eigenvalue's distribution  
    # compared to the Marchenko–Pastur distribution of random matrices; 
    # differences suggest biological significance of distances
    
    return eigenvalues

def largest_connected_component(L): # returns largest cc of laplacian matrix

    n_components, labels = csgraph.connected_components(L)

    if n_components == 1:
        return L

    # find largest component
    counts = np.bincount(labels)
    largest_label = np.argmax(counts)

    keep = labels == largest_label
    return L[keep][:, keep]

if __name__ == '__main__':
    main()

