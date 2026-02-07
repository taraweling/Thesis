import csv # do i need this?
from pyvis.network import Network as net
import networkx as nx
import matplotlib.pyplot as plt


def kmeans(adjlist): # UNFINISHED # struggles with clusters of different densities (aka evens out size)
    
    
    
    return

def calculate_clustering_coeff(adjlist): # input = adjlist like the one from A
    # output = dict of (node:str, clustering_coeff:float)
    C = {}
    clustering_coeff = 0
        
    for i in adjlist:
        
        # find total number of edges by adding all values of the input dictionary for denominator
        degree = len(adjlist[i])
        emax = degree * (degree-1)
        
        # check if denominator = 0 and equation is undefined
        if emax == 0:
            C[i] = 0.0 # because float
            continue # skip to next iteration of loop
        
        # numerator = check if neighbors of v are connected to each other
        
        edge = 0 # tracker
        # check every possible value of pair (u,w)
        for u in adjlist[i]:
            for w in adjlist[i]:
                #if w in adjlist[u]: # could also check if u in w's adj list for directed graphs         
                if u in adjlist and w in adjlist and w in adjlist[u]: # chat gpt corrected my original "if w in adjlist[u]"
                    edge += 1
        
        # formula 
        clustering_coeff = edge/emax # shouldn't it be edge/2 to correct for double counting?
        
        # C[i] is key
        C[i] = clustering_coeff
    
    return C

# I received help from CS drop-in tutor Doran Penner for the below fn
def calculate_closeness(adjlist): #input adjlist from A, output dict of (node, closeness) pairs
    S = {} # dict to return
    closeness = 0.0
    #  closeness centrality = reciprocal of the sum of all nodes' shortest paths to v.
    # also see: len of adjlist - 1 / sum of shortest paths to each nodes

    # initiate and populate distances dictionary
    
    nodelist = [] # all nodes excluding i 
    
    for i in adjlist:
        # 
        for n in adjlist: # make nodelist
            if n != i:
                nodelist.append(n)
    
        # length of shortest path from u to v = apply BFS          
        shortestpathdict = {} # shortest paths dictionary
        for k in nodelist:
            shortestpathdict[k] = float('inf')
            
        shortestpathdict[i] = 0 # set starting node to 0 
        
        # edited Anna's BFS code from lab2_sol.py to return shortest paths from s to all other nodes
        Q = [i] # initialize and populate queue of nodes to explore
        
        while len(Q) != 0:
            exploring = Q.pop(0)
            for neighbor in adjlist[exploring]: ## for each of exploring's neighbors...
                if shortestpathdict[neighbor] == float('inf'): # update if we haven't reached that neighbor yet
                    shortestpathdict[neighbor] = shortestpathdict[exploring] + 1
                    Q.append(neighbor)
        
        # hi: shortestpathdict bw non-i and everything else
        ## calc closeness of i

        numerator = (len(adjlist))-1.0
        denominator = 0.0
        # with dictionary of distances, sum up all values from 'u's to v 
        denominator = sum(shortestpathdict.values())
    
        # apply formula
        closeness = numerator/denominator
        
        S[i] = closeness

    return S