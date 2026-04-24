import statistics #obvious package is obvious
import math #obvious package #2 is obvious
from collections import deque
# by thomas

# Functions will be standardized as documented in graph_algos.py. 
## General order of inputs is adjlist, adjmatrix, nodelist, edgelist (will be edited if I change my mind), target node
# data_summary_statistics returns list mean, median, mode, min, max

def unweight_adjlist(weighted_adjlist, nodelist):
    adjlist = {}
    for node in nodelist:
        adjlist[node] = []
        for neighbor in weighted_adjlist[node]:
                adjlist[node].append(neighbor[0])
    return adjlist
#just chop up unweight_adjlist as you see fit

def node_degree_stats(adjlist):
    degree_list = []
    for node in adjlist.keys():
        degree_list += len(adjlist[node])
    return data_summary_statistics(degree_list) 
#returns mean, median, mode, min, max in order

def distance_stats(weighted_adjlist): # THESE TWO TAKE THE SAME INPUT
    return data_summary_statistics(get_edge_weights(weighted_adjlist))

def weight_stats(weighted_adjlist): # PROVIDE DATA AS GIVEN BY MILES' ALGO
    return data_summary_statistics(get_edge_weights(reciprocals(weighted_adjlist)))

def weighted_node_degree(weighted_adjlist, node): #these work with either kind of weighting
    node_degree = 0
    for neighbor, weight in weighted_adjlist[node]:
        node_degree += weight
    return node_degree

def weighted_node_degree_stats(weighted_adjlist): #these work with either kind of weighting
    degree_list = []
    for node in weighted_adjlist.keys():
        degree_list.append(weighted_node_degree(weighted_adjlist, node))
    return data_summary_statistics(degree_list) 

def eccentricity_centrality_distribution(weighted_adj_list): 
    eccentricity_list = []
    for node in weighted_adj_list.keys():
        D = single_dijkstra(weighted_adj_list, node)
        eccentricity_list.append(1/max(D.values()))
    return data_summary_statistics(eccentricity_list)

def eccentricity_centralities(weighted_adj_list):
    eccentricity_dict = {}
    for node in weighted_adj_list.keys():
        D = single_dijkstra(weighted_adj_list, node)
        vals = []
        for d in D.values():
            if d != float('inf') and d != 0:
                vals.append(d)

        if len(vals) == 0:
            eccentricity_dict[node] = 0.0
        else:
            eccentricity_dict[node] = 1 / max(vals)
            
        """vals = D.values()
        eccentricity_dict[node] = float(1/(max(vals)))"""
    return eccentricity_dict

def closeness_centrality_distribution(weighted_adj_list): 
    closeness_list = []
    for node in weighted_adj_list.keys():
        D = single_dijkstra(weighted_adj_list, node)
        closeness_list.append(1/(sum(D.values())))
    return data_summary_statistics(closeness_list)

def closeness_centralities(weighted_adj_list): 
    closeness_dict = {}
    for node in weighted_adj_list.keys():
        D = single_dijkstra(weighted_adj_list, node)
        vals = []
        # below = tara's addition to ignore inf
        for d in D.values():
            if d != float('inf') and d != 0:
                vals.append(d)

        if len(vals) == 0:
            closeness_dict[node] = 0.0
        else:
            closeness_dict[node] = 1 / sum(vals)
    return closeness_dict

def betweenness_centrality_stats(weighted_adj_list, normalized):
    counts = {node:0 for node in weighted_adj_list.keys()}
    for node in weighted_adj_list.keys():
        betweenness_updater(weighted_adj_list, node, counts)
    betweenness_values = counts.values()
    betweenness_values = list(betweenness_values)
    for i in range(len(betweenness_values)):
        betweenness_values[i] = betweenness_values[i]/(len(weighted_adj_list.keys())*(len(weighted_adj_list.keys())-1))
    maxi = max(betweenness_values)
    mini = min(betweenness_values)
    if normalized == True:
        normalized_betweenness = []
        for i in range(len(betweenness_values)):
            normalized_betweenness.append((betweenness_values[i] - mini)/(maxi - mini))
            return data_summary_statistics(normalized_betweenness)
    else:
        return data_summary_statistics(betweenness_values)

def betweenness_centralities(weighted_adj_list, normalized):
    counts = {}
    for node in weighted_adj_list.keys():
        counts[node] = 0
    for node in weighted_adj_list.keys():
        betweenness_updater(weighted_adj_list, node, counts)
    if normalized == True:
        maxi = max(counts.values())
        mini = min(counts.values())
        for i in counts.keys():
           counts[i] = ((counts[i] - mini)/(maxi - mini))
        return counts
    else:
        return counts

def calculate_clustering_coeff(adjlist): 
    C = {}
    for node, neighbors in adjlist.items():
        if len(neighbors) < 2:
            C[node] = 0
            continue
        neighbor_edge_count = 0
        for u in neighbors:
            for v in neighbors:
                neighbor_edge_count += 1 if u in adjlist[v] else 0
        C[node] = neighbor_edge_count/(len(neighbors)*(len(neighbors)-1))
    return C

def graph_cluster_coeff(adjlist):
    nodelist = adjlist.keys()
    adjmat = adjlist_to_adjmat(adjlist)
    triangles_made = 0
    for i in nodelist:
        for j in nodelist:
            for k in nodelist:
                if i<j and j<k:
                    triangles_made += (adjmat[i][j]*adjmat[j][k]*adjmat[i][k])
    possible_triangles = math.comb(len(nodelist), 3)
    graph_coeff = triangles_made/possible_triangles
    return graph_coeff

# helper functions placed down here for readability
def single_dijkstra_pred(weighted_adj_list, s): # make sure not to double count edges
    _, pred, _ = bellman_ford(weighted_adj_list, s)
    return pred

def betweenness_updater(weighted_adj_list, start, counts):
    pred = single_dijkstra_pred(weighted_adj_list, start)
    for node in pred.keys():
        current_node = node
        while current_node is not None:
            counts[current_node] +=1
            current_node = pred[current_node]

def single_dijkstra(weighted_adj_list, s): # make sure not to double count edges
    # keep this function for backwards compatibility; now delegated to Bellman-Ford
    D, _, _ = bellman_ford(weighted_adj_list, s)
    return D

def bellman_ford(weighted_adj_list, s):
    """
    Bellman-Ford on a weighted adjacency list.
    Returns:
      dist: shortest-distance estimates from s
      pred: predecessor map for one shortest-path tree
      neg_cycle_nodes: nodes on or reachable from a negative cycle reachable from s
    """
    nodes = set(weighted_adj_list.keys())
    edges = []

    for u, neighbors in weighted_adj_list.items():
        for v, w in neighbors:
            nodes.add(v)
            edges.append((u, v, float(w)))

    dist = {node: float('inf') for node in nodes}
    pred = {node: None for node in nodes}
    if s not in dist:
        dist[s] = 0.0
        pred[s] = None
    else:
        dist[s] = 0.0

    # Relax edges |V|-1 times.
    for _ in range(max(0, len(nodes) - 1)):
        updated = False
        for u, v, w in edges:
            if dist[u] == float('inf'):
                continue
            candidate = dist[u] + w
            if candidate < dist[v]:
                dist[v] = candidate
                pred[v] = u
                updated = True
        if not updated:
            break

    # Detect negative cycles reachable from s.
    neg_cycle_nodes = set()
    for u, v, w in edges:
        if dist[u] == float('inf'):
            continue
        if dist[u] + w < dist[v]:
            neg_cycle_nodes.add(u)
            neg_cycle_nodes.add(v)

    # Propagate negative-cycle effect to all downstream nodes.
    if neg_cycle_nodes:
        queue = deque(neg_cycle_nodes)
        while queue:
            curr = queue.popleft()
            for neighbor, _ in weighted_adj_list.get(curr, []):
                if neighbor not in neg_cycle_nodes:
                    neg_cycle_nodes.add(neighbor)
                    queue.append(neighbor)

        for node in neg_cycle_nodes:
            dist[node] = -float('inf')
            if node != s:
                pred[node] = None

    return dist, pred, neg_cycle_nodes

def reciprocals(weighted_adjlist): #swaps weighting between distance and Weight so keep track in main
    reciprocal_weights = {}
    for node in weighted_adjlist.keys():
        reciprocal_weights[node] = []
        for neighbor, weight in weighted_adjlist[node]:
            reciprocal_weights[node].append([neighbor, max(get_edge_weights(weighted_adjlist)) - weight + 1]) 
    return reciprocal_weights

def data_summary_statistics(list): #returns mean, median, mode, min, max of an unsorted list
    list.sort() #literally just for peace of mind
    mean = statistics.mean(list)
    median = statistics.median(list)
    mode = statistics.mode(list)
    minimum = min(list)
    maximum = max(list)

    return list, mean, median, mode, minimum, maximum #in order

def adjmat_to_adjlist(adj_mat, nodelist): #self explanatory
    adjlist = {}
    for i in range(len(nodelist)):
        adjlist[nodelist[i]] = []
        for j in range(len(nodelist)):
            if adj_mat[i][j]==1:
                adjlist[nodelist[i]].append(nodelist[j])
    return adjlist

# ez weighted algo ez 
def read_edges_weighted(edgefile):
    adjlist = {} 
    with open(edgefile, "r", encoding="utf-8") as f: ## short for file
        lines = f.readlines() ## list of lines
        for line in lines:
            key, value, weight = line.strip().split("	") ## strip and split my beloveds
            if key not in adjlist: ## if key isn't in the dictionary, make an empty set as the value before we append
                adjlist[key] = []
            adjlist[key].append([value, weight])
            if value not in adjlist: ## identical code to above comment
                adjlist[value] = []
            adjlist[value].append([key, weight])
    return adjlist

def get_nodes(adjlist): # returns nodelist
    nodelist = list(adjlist)
    return nodelist

def get_edge_weights(weighted_adjlist): #self explanatory
    edge_weights = []
    for node in weighted_adjlist.keys():
        for neighbor, weight in weighted_adjlist[node]:
            if neighbor < node:
                #print("edge weight thomas",weight,int(weight)) 
                edge_weights.append(weight) # changed from += to .append()
                
    return edge_weights

# obvious functions are obvious, shamelessly stolen from p2 and likely obsolete ty miles for making it output an adjlist
def read_edges(edgefile):
    adjlist = {} 
    with open(edgefile, "r", encoding="utf-8") as f: ## short for file
        lines = f.readlines() ## list of lines
        for line in lines:
            key, value = line.strip().split("   ") ## strip and split my beloveds
            if key not in adjlist: ## if key isn't in the dictionary, make an empty set as the value before we append
                adjlist[key] = []
            adjlist[key].append(value)
            if value not in adjlist: ## identical code to above comment
                adjlist[value] = []
            adjlist[value].append(key)
    nodelist = adjlist.keys()
    return nodelist, adjlist
# read_edges returns both a nodelist AND an adjlist

"""def adjlist_to_adjmat(nodelist, adjlist):
    adjmat = [[0 for node2 in range(len(nodelist))] for node1 in range(len(nodelist))]
    # we define a matrix of all zeroes
    for node in nodelist:
        for neighbor in adjlist[node]:
            adjmat[node][neighbor] = 1
            adjmat[neighbor][node] = 1
    return adjmat"""

def adjlist_to_adjmat(adjlist): # redid above fn to input weighted adjlist using chatgpt
    # n = len(adjlist)
    # #adjmat = [[0] * n for _ in range(n)]
    # adjmat = [[0.0 for _ in range(n)] for _ in range(n)]
    
    # for node, neighbors in adjlist.items():
    #     for neighbor, weight in neighbors:
    #         adjmat[node][neighbor] = weight
    #         adjmat[neighbor][node] = weight 
    #     # for neighbor in neighbors:
    #     #     print(node,"+",neighbor)
    #     #     adjmat[node][neighbor[0]] = 1

    # collect all nodes (keys + neighbors)
    nodes = set(adjlist.keys())
    for neighbors in adjlist.values():
        for neighbor, _ in neighbors:
            nodes.add(neighbor)

    nodes = sorted(nodes)
    node_to_idx = {node: i for i, node in enumerate(nodes)}

    n = len(nodes)

    # initialize adjacency matrix
    adjmat = [[math.inf for _ in range(n)] for _ in range(n)]
    for i in range(n):
        adjmat[i][i] = 0.0

    # fill matrix
    for node, neighbors in adjlist.items():
        i = node_to_idx[node]
        for neighbor, weight in neighbors:
            j = node_to_idx[neighbor]
            adjmat[i][j] = weight
            adjmat[j][i] = weight

    return adjmat


def normalize_to_unit_values(values): #normalize to 0-1 linearly
    maxi = max(values)
    normalized = []
    for val in values:
        normalized.append(val/maxi)
    return normalized
        
def search_array(array, value):
    for n in range(len(array)):
        if array[n] == value:
            return n
        
