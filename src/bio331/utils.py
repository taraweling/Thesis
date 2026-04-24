import sys
import copy
import random
import math
from pyvis.network import Network

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

def read_foodweb_nodes(foodfile):
    nodes = {}
    
    with open(foodfile) as fin:
        
        for line in fin:
            
            row = line.strip().split()
            
            nodes[row[0]] = [row[1]] # taxon name
            nodes[row[0]].append(row[2]) # classification
    return nodes

def read_foodweb_edges(infile): # .edited from my code for p2
    edgelist = {} # return dictionary of sets
    
    # below three lines taken from lab 1 solutions
    unique = [] # unique node names tracker
    edges = [] # list of edges
    
    with open(infile) as fin: 
        
        for line in fin:
            
            row = line.strip().split()
                
            for node in row: # add keys to dictionary, then search through edgelist for values
                # makes a list of unique nodes
                if str(node) not in unique:
                    
                    edgelist[str(node)] = {}
                    
                    unique.append(str(node))
    
    # for key in adjlist, if key matches with node then add as value to said key
        for u, v in edges:
            
            if u in edgelist.keys() and edgelist[u]: #and adjlist[u]: # updates 
                edgelist[u].update([v]) # make list to not separate each char
                # because update takes iterator
                
            # if there's already a value for the u key, create nested dict value instead of overwriting
            elif u in edgelist.keys():
                edgelist[u] = {v}
            
            if v in edgelist.keys() and edgelist[v]:
                edgelist[v].update([u])
                
            # elif there's already a value for the u key, add to nested dict value instead of overwriting
            elif v in edgelist.keys():
                edgelist[v] = {u} 
    return edgelist

# incomplete
def viz_foodweb(nodes,edges,partition,outfile):
    """
    Visualize a food web and write it to an HTML file.
    :param: nodes - list or set of nodes
    :param: edges - list of 2-elemenåt lists.
    :param: nodenames -
    :param: partition - list of lists representing a partition of the nodes to color.
    :param: outfile - string outfile that ends in '.html'
    :returns: None
    """
    #

    """    # get colors for partitions
    node_colors = {}
    for cluster in partition:
        color = rgb_to_hex(random.random(),random.random(),random.random())
        for n in cluster:
            node_colors[n] = color
    """
    
         
    """
    for n in nodes: # add nodes
        if 'Benthic' in nodenames:
            G.add_node(n,label=nodenames[0].values(),shape='square')
        if 'Pelagic' in nodenames:
            G.add_node(n,label=nodenames[0].values(),shape='triangle')
        if 'Unknown' in nodenames:
            G.add_node(n,label=nodenames[0].values(),shape='diamond')
    """
    
    # create graph
    # annotate node label according to taxon name
    G = Network() 
    node_shapes = {}
 
    # get shapes for partitions
    node_shapes = {}
    for name in nodes:
        shape = ('square','triangle','diamond')
        for n in name:
            node_shapes[n] = shape
            
    # get labels for partitions
    node_labels = {}
    print(i[1] for i in nodes)
    for name in nodes:
        print(name)
        for n in name:
            node_labels[n] = n

    G = Network() # create graph
    for n in nodes: # add nodes
        if 'Benthic' in nodes[n]:
            G.add_node(n,label=n,shape='square')
        if 'Pelagic' in nodes[n]:
            G.add_node(n,label=n,shape='triangle')
        if 'Unknown' in nodes[n]:
            G.add_node(n,label=n,shape='diamond')
    
    
    for u in edges:
        for v in edges[u]:
            G.add_edge(u,v)

    G.toggle_physics(True) 
    G.show_buttons(filter_=['physics'])
    # taken from lab 2, 

    G.write_html(outfile)
    print('wrote to', outfile)
    return
