import utils
import random
import sys
import math
from pyvis.network import Network

def main():
    ## Task A: Read the Edge File
    nodes, edgelist, adjlist = read_edges('files/example-edges.txt')
    
    print('nodes:', nodes)
    print('edges:', edgelist)
    print('adjlist', adjlist)
    
    ## Task B: Edge Betweenness

    ### B1. Calculate the edge betweenness for one edge.
    paths = all_pairs_shortest_paths(adjlist) 
    print(paths,'\n')
    print(single_edge_betweenness(nodes,paths,'v3','v5')) # should be 15.0
    print(single_edge_betweenness(nodes,paths,'v5','v7')) # should be 7.5

    ### B2. Calculate edge betweenness for all edges in the graph.
    B = edge_betweenness(nodes,edgelist,paths)
    for key in B:
        print(key,B[key])

    ## Task C: Implement the Girvan-Newman Algorithm
    partitions = GN(nodes,edgelist,adjlist)
    for p in partitions:
        print(p)
     
    """
    testedge = [['a','b'],['b','c'],['d','a']]
    utils.remove_from_edgelist(['a','b'],testedge)
    print(testedge)
    utils.remove_from_edgelist(['a','d'],testedge)
    print('TESTING:',testedge)
    """
    
    ## Task D: Post the Example Graph Colored by Communities
    nodes, edgelist, adjlist = read_edges('files/example-edges.txt')
    utils.viz_example(nodes,edgelist,partitions[2],'3groups.html') # I'm not really sure what to do beyond this
    utils.viz_example(nodes,edgelist,partitions[3],'4groups.html')
    utils.viz_example(nodes,edgelist,partitions[3],'5groups.html')
    
    ## Task E: Apply and Visualize Girvan-Newman on the Chesapeake Bay Food Web
    ### read in foodweb from files/foodweb-edges.txt and foodweb-nodes.txt 
    food_nodes, food_edgelist, food_adjlist = read_edges('files/foodweb-edges.txt')
    print(food_nodes,'\n',food_edgelist,'\n',food_adjlist)
    
    # foodweb-nodes file is where taxa (numbers) have ID
    node_names = utils.read_nodes('files/foodweb-nodes.txt')
    
    ### run GN to generate all partitions
    food_partitions = GN(food_nodes,food_edgelist,food_adjlist)
    for p in food_partitions:
        print(p)
        
    ### visualize chesapeake bay food web = incomplete, calls are commented out
    foodweb_edgelist = utils.read_foodweb_edges('files/foodweb-edges.txt')
    foodweb_nodes = utils.read_foodweb_nodes('files/foodweb-nodes.txt')

    utils.viz_foodweb(foodweb_nodes,foodweb_edgelist,partitions[3],'foodweb-2clusters.html') 
    #utils.viz_foodweb(nodes,edgelist,node_names,partitions[4],'foodweb-3clusters.html')
     
    ## Task F: Compare the Unweighted Clustering to the Paper Result
    """
    Similarities and differences compared to public results:
    """
    
    return # done with main()


def read_edges(infile): # .edited from my code for p2
    nodes = set() # convert to list later
    edgelist = [] # 2 element lists
    
    adjlist = {} # adjlist dict to return
    
    with open(infile) as fin:
        for line in fin:
            u, v = line.strip().split()

            edgelist.append((u, v))
            nodes.update([u, v])

            if u not in adjlist:
                adjlist[u] = set()
                
            if v not in adjlist:
                adjlist[v] = set()

            adjlist[u].add(v)
            adjlist[v].add(u)

    return list(nodes), edgelist, adjlist

## taken from lab 5 solutions, didn't add to utils because it broke test_run.py
def all_pairs_shortest_paths(adjlist):
    nodes = adjlist.keys() # all nodes in the graph
    paths = {}
    for n in nodes:
        paths[n] = {} # dictionary of dictionaries for every node.

        D,pred = utils.bellman_ford_ties(adjlist,n) # Bellman-Ford predecessor ties
        

        for node in D.keys():
            if node != n and math.isfinite(D[node]): # skips source and unreachable/neg-cycle nodes
                paths[n][node] = utils.get_paths(pred,node,[[node]])

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

# Task B2 dict {edges=tuples : edge betweenness centrality=floats}
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
        utils.remove_from_adjlist(forward,adjlist)
        utils.remove_from_edgelist(forward,edgelist)
        
        utils.remove_from_adjlist(reverse,adjlist)
        # utils.remove_from_edgelist(reverse,edgelist)
        print(adjlist,edgelist)

        # was not getting (v3, v5) but (v1,v2) = 1 
        ## copied my code into chatgpt and I had forgotten to check for reverse (v,u) at line 125

        ### C4. If removing an edge divided a group, make a new partition.
        ## call conncomp = returns sorted list of nodes in the same connected component as u
        print('1',forward[0],'\n','2',forward[1])
        list1 = utils.conn_comp(adjlist,forward[0]) # list of nodes in same connected component as v3
        list2 = utils.conn_comp(adjlist,forward[1]) # list of nodes in same connected component as v5

        if list1 != list2: # new partition
            new = utils.split_partition(partitions[-1],list1,list2)
            partitions.append(new)
        
        ### C5. Repeat from Step 2 until no edges remain to remove.
        edge_between  = edge_betweenness(nodes,edgelist,all_pairs_shortest_paths(adjlist))

    return partitions

# keep this at the bottom of the file.
if __name__ == '__main__':
    main()
