import sys
import copy
import random
import math
from pyvis.network import Network

## taken from lab 5 solutions, didn't add to utils because it broke test_run.py
def all_pairs_shortest_paths(adjlist):
    nodes = adjlist.keys() # all nodes in the graph
    paths = {}
    for n in nodes:
        paths[n] = {} # dictionary of dictionaries for every node.

        D,pred = bellman_ford_ties(adjlist,n) # Bellman-Ford predecessor ties
        

        for node in D.keys():
            if node != n and math.isfinite(D[node]): # skips source and unreachable/neg-cycle nodes
                paths[n][node] = get_paths(pred,node,[[node]])

    return paths

# Girvan-Newman betweenness problem
## adapted from lab 5 solutions, added to main instead of utils because of test_run.py's specifications
def single_edge_betweenness(nodes,paths,u,v):
    B = 0
    
    for i in range(len(nodes)):
        for j in range(i+1, len(nodes)):
            s = nodes[i]
            t = nodes[j]
            if t in paths[s]:
                numerator = 0
                for p in paths[s][t]:
                    for k in range(len(p)-1):
                        if (p[k] == u and p[k+1] == v) or (p[k] == v and p[k+1] == u):
                            numerator += 1
                B += numerator / len(paths[s][t])

    return B

#{edges=tuples : edge betweenness centrality=floats}
def edge_betweenness(nodes, edges, paths): # paths calls single_edge_betweenness
    B = {}
    for i in edges: # call single_edge for every single edge
        B[i] = single_edge_betweenness(nodes,paths,i[0],i[1])
    #x = single_edge_betweenness()
    return B

def GN(nodes,edgelist,adjlist): # outputs list of partitions where each is a list of node groups (lists)
    ### C1. Assign all nodes to a single cluster.
    partitions = [[nodes]]
    while len(edgelist) > 0:
        ### C2. Calculate the _edge betweenness_ of all edges in the network.  
        paths = all_pairs_shortest_paths(adjlist)
        edge_bw = edge_betweenness(nodes,edgelist,paths)
        
        ### C3. Remove the edge with the highest betweenness
        highest = max(edge_bw, key = edge_bw.get) # googled ways to index dictionaries
        forward = list(highest)
        # reverse
        reverse = list(highest)
        reverse.reverse()
        
        ## remove forward  edge from graph
        remove_from_adjlist(forward,adjlist)
        remove_from_edgelist(forward,edgelist)
        
        remove_from_adjlist(reverse,adjlist)
        # utils.remove_from_edgelist(reverse,edgelist)
        print(adjlist,edgelist)

        # was not getting (v3, v5) but (v1,v2) = 1 
        ## copied my code into chatgpt and I had forgotten to check for reverse (v,u) at line 125

        ### C4. If removing an edge divided a group, make a new partition.
        ## call conncomp = returns sorted list of nodes in the same connected component as u
        print('1',forward[0],'\n','2',forward[1])
        list1 = conn_comp(adjlist,forward[0]) # list of nodes in same connected component as v3
        list2 = conn_comp(adjlist,forward[1]) # list of nodes in same connected component as v5

        if list1 != list2: # new partition
            new = split_partition(partitions[-1],list1,list2)
            partitions.append(new)
        
        ### C5. Repeat from Step 2 until no edges remain to remove.
        edge_between  = edge_betweenness(nodes,edgelist,all_pairs_shortest_paths(adjlist))

    return partitions


#######
## Recursive function that calculates all possible shortest paths
## for a node based on the predecessors from shortest_paths() function.
## It is called for the LAST node in the path and traces paths backwards.
## Inputs: predecessors dictionary
## Inputs: current node n to process
## Inputs: current paths list (paths from the end traced back to n)
## Returns: List of paths, where each path is a list of nodes
#######
def get_paths(predecessors,n,curr_paths):
    #print(n,curr_paths)
    # base case
    # we've reached the first node; we're done.
    if predecessors[n] == None:
        return curr_paths

    # there's still at least one node to backtrack.
    # for each predecessor, prepend the predecessor and call
    # get_paths() for the predecessor.
    new_paths = []
    for pred in predecessors[n]: 
        these_paths = []
        for i in range(len(curr_paths)):
            these_paths.append([pred]+curr_paths[i])
        new_paths+=get_paths(predecessors,pred,these_paths)
    return new_paths

#######
## Removes an edge from an edgelist.
#######
def remove_from_edgelist(edge,edge_list):
    # get edge from edge list. Could be (u,v) or (v,u)
    u = edge[0]
    v = edge[1]
    options = [(u,v),(v,u),[u,v],[v,u]]
    # 
    for e in options:
        if e in edge_list:
            edge_list.remove(e)
    ## returns nothing - removes in-place
    return

#######
## Removes an edge from an adjacency list.
#######
def remove_from_adjlist(edge,adj_list):
    u = edge[0]
    v = edge[1]
    adj_list[u].remove(v)
    ## returns nothing - removes in-place
    return

#######
## Given a partition, a current group, and the group divided in two,
## returns a new partition with the current group removed and
## the two split groups added.
#######
def split_partition(prev_partition,split1,split2):
    # if statement to check if partition is not a list of lists
    return [split1, split2]
    """if any(isinstance(item, list) for item in prev_partition): # line generated by asking chatgpt to write a function to check if a list is a list of lists
        partition = copy.deepcopy(prev_partition) # make a copy of partition
        list_to_remove = []
        before = split1+split2
        for p in partition:
            if sorted(p)==sorted(before):
                list_to_remove = p
    
        partition.remove(list_to_remove)
        partition.append(sorted(split1))
        partition.append(sorted(split2))
        return partition
    else:
        return [split1, split2]"""

#######
## Given an adjacency list and a node u, 
## returns a sorted list of nodes that are 
## in the same connected component as u.
#######
def conn_comp(adjlist,u):
    # initialize a queue Q and a set of seen nodes
    Q = [u]
    seen = set()
    seen.add(u)

    while len(Q) > 0: # while there's still a node to explore...
        exploring = Q.pop(0) # remove the FIRST node from Q
        for neighbor in adjlist[exploring]:
            if neighbor not in seen: # unexplored
                seen.add(neighbor) # add the neighbor to the seen node set
                Q.append(neighbor) # append the neighbor to Q
    return sorted(list(seen))


""" Below four functions taken from utils.py for lab 5 """
#######
## Computes the shortest paths from n to all other
## nodes, also returns predecessors of nodes, keeping ties.
## Inputs: adj_list (dictionary) - adjacency list
## Inputs: n (string) - node
## Returns: distance dictionary (u:distance) for every node u
## Returns: predecessors dictionary (u:[pred list]) for every node u
#######
def bellman_ford_ties(adj_list, n):
    """
    Bellman-Ford with predecessor tie tracking.
    Supports unweighted adjacency lists (neighbors only) and weighted ones ([neighbor, weight]).
    """
    nodes = set(adj_list.keys())
    edges = []

    for u, neighbors in adj_list.items():
        for item in neighbors:
            if isinstance(item, (list, tuple)) and len(item) >= 2:
                v = item[0]
                w = float(item[1])
            else:
                v = item
                w = 1.0
            nodes.add(v)
            edges.append((u, v, w))

    dist = {node: float('inf') for node in nodes}
    predecessors = {node: [] for node in nodes}
    dist[n] = 0.0
    predecessors[n] = None

    for _ in range(max(0, len(nodes) - 1)):
        updated = False
        for u, v, w in edges:
            if dist[u] == float('inf'):
                continue
            candidate = dist[u] + w
            if candidate < dist[v]:
                dist[v] = candidate
                predecessors[v] = [u]
                updated = True
            elif math.isclose(candidate, dist[v]) and predecessors[v] is not None and u not in predecessors[v]:
                predecessors[v].append(u)
        if not updated:
            break

    neg_cycle_nodes = set()
    for u, v, w in edges:
        if dist[u] != float('inf') and dist[u] + w < dist[v]:
            neg_cycle_nodes.add(v)
            neg_cycle_nodes.add(u)

    if neg_cycle_nodes:
        frontier = list(neg_cycle_nodes)
        while frontier:
            curr = frontier.pop()
            for item in adj_list.get(curr, []):
                if isinstance(item, (list, tuple)) and len(item) >= 1:
                    nxt = item[0]
                else:
                    nxt = item
                if nxt not in neg_cycle_nodes:
                    neg_cycle_nodes.add(nxt)
                    frontier.append(nxt)
        for node in neg_cycle_nodes:
            dist[node] = -float('inf')
            if node != n:
                predecessors[node] = []

    return dist, predecessors

def shortest_paths_ties(adj_list, n):
    return bellman_ford_ties(adj_list, n)


## RGB to Hex function - copied from Lab 2
def rgb_to_hex(red,green,blue): # pass in three values between 0 and 1
  maxHexValue= 255  ## max two-digit hex value (0-indexed)
  r = int(red*maxHexValue)    ## rescale red
  g = int(green*maxHexValue)  ## rescale green
  b = int(blue*maxHexValue)   ## rescale blue
  RR = format(r,'02x') ## two-digit hex representation
  GG = format(g,'02x') ## two-digit hex representation
  BB = format(b,'02x') ## two-digit hex representation
  return '#'+RR+GG+BB

## visualize the graph.
def viz_example(nodes,edges,partition,outfile):
    """
    Visualize a graph and write it to an HTML file.
    :param: nodes - list or set of nodes
    :param: edges - list of 2-element lists.
    :param: partition - list of lists representing a partition of the nodes to color.
    :param: outfile - string outfile that ends in '.html'
    :returns: None
    """
    # Refer to Lab 1 for instructions about visualizing a graph.

    # get colors for partitions
    node_colors = {}
    for cluster in partition:
        color = rgb_to_hex(random.random(),random.random(),random.random())
        for n in cluster:
            node_colors[n] = color

    G = Network() # create graph
    for n in nodes: # add nodes
        G.add_node(n,label=n,color=node_colors[n])
    for u,v in edges: # add edges
        G.add_edge(u,v) 

    G.toggle_physics(True) 
    G.show_buttons(filter_=['physics'])

    G.write_html(outfile)
    print('Saved file as',outfile)

    return


## read nodes from p2
def read_nodes(nodefile): # input = file, output = dict of (node, label) pairs for every node
    nodes = {}

    # skip header file
    with open(nodefile) as fin: 
        
        for line in fin:
            
            row = line.strip().split()
            
            if row[0] != 'Node': # ignores header
                nodes[row[0]] = row[1]
                        
    return nodes

def viz_foodweb(nodes,nodenames,edges,partition,outfile):
    """
    Visualize a food web and write it to an HTML file.
    :param: nodes - list or set of nodes
    :param: edges - list of 2-element lists.
    :param: partition - list of lists representing a partition of the nodes to color.
    :param: outfile - string outfile that ends in '.html'
    :returns: None
    """
    ## annotate node label according to taxon name

    # get colors for partitions
    node_colors = {}
    for cluster in partition:
        color = rgb_to_hex(random.random(),random.random(),random.random())
        for n in cluster:
            node_colors[n] = color

    G = Network() # create graph
    for n in nodes: # add nodes
        G.add_node(n,label=n,color=node_colors[n])
    for u,v in edges: # add edges
        G.add_edge(u,v) 

    G.toggle_physics(True) 
    G.show_buttons(filter_=['physics'])

    G.write_html(outfile)
    print('Saved file as',outfile)

    return
