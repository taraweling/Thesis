import csv # do i need this?
from pyvis.network import Network as net
import networkx as nx
import matplotlib.pyplot as plt

# use below fn on deg_grn_tfsonly: 
def regulator_detection(grn, disorderlist):

    """
    Identify TFs whose high-weight regulatory edges target genes with large expression shifts = regulator drivers in PANDA network.

    Driver score per TF: sum(weight * abs(log2fc_target))
    = estimates how strongly a transcription factor's regulatory edges align with observed expression changes.

    Output:
    {'TF':..., 'targets':..., 'deg_targets':..., 'driver_score':..., 'mean_target_log2fc':...}
    """
    
    # Build gene: disorder metadata dictionary
    gene_to_disorders = {}

    for rec in disorderlist:
        gene_id = rec[0]
        info = tuple(rec[1:])
        gene_to_disorders.setdefault(gene_id, []).append(info)

    degset = set(gene_to_disorders.keys())

    tf_results = []

    for tf, targets in grn.items():

        score = 0
        deg_targets = 0
        total_targets = 0
        log2fcs = []

        for gene, weight in targets:

            total_targets += 1

            if gene not in degset:
                continue

            for disorder_info in gene_to_disorders[gene]:
                log2fc = disorder_info[4]
                score += weight * abs(log2fc)
                log2fcs.append(log2fc)
                deg_targets += 1

        if deg_targets == 0:
            continue

        mean_fc = sum(log2fcs) / len(log2fcs)
        
        normalized_score = score / total_targets 
        
        tf_results.append({
            "TF": tf,
            "targets": total_targets,
            "deg_targets": deg_targets,
            "driver_score": score, # = total regulatory influence
            "normalized_score": normalized_score, # # = average influence per target
            "mean_target_log2fc": mean_fc})
    
        tf_results.sort(key=lambda x: x["driver_score"], reverse=True)

    return tf_results


def edgeweight_summary(adjlist,*condition:None): #works for any grn 
    """
    Input:
        Any adjacency list structured as:
            {TF: (gene1, regweight1), (gene2, regweight2)]}
    
    Output:
        printable dict containing average edge weight counts for:
            - DETFs (TF nodes)
            - DEGs (targets)
            - DEXs (both tfs and gene targets)
            
    """
    tf_avg = [] # tf avgs
    g_avg = []
    degs_avg = [] # list of all DEG weights 
    detfs_avg = []
    
    for tf, edges in adjlist.items():
        
        for edge in edges: # loops through each gene
            weight = edge[1]
            if len(edges) > 2:
                degs_avg.append(weight)
            else:
                g_avg.append(weight) # 
                
        tf_avg.append(sum(degs_avg) / len(degs_avg)) # adds average deg per tf
        
        if 
    summary = { 
        "DEGs_average":  sum(degs_avg) / len(degs_avg), # averages all degs
        "TFs_avg": sum(tf_avg) / len(tf_avg), # average TF
        "DETFs_avg": sum(),
        "G_avg": sum(g_avg) / len(g_avg)
        }
    
    print(summary) # remove later 
    return summary

# what about a function to look at differential TFs and the genes they regulate? 
# would require 

def log2fc_summary(adjlist): # only applies to degs = whatever has 7 values 
    #per tf (which does not have to be differentially regulated)
    
    """
    Input:
        DEG-GRN adjacency list structured as :
            {TF: (DEG1, etc), (DEG2, etc)]}
            
        Figure out later if I should cut TFs down to DETF...
            {TF: (geneID, edgeweight, disorder, study, year, tissue, log2fc, pval)}

    Output:
        dict containing sign counts for:
            - DEGs (targets)
            - DETFs (TF nodes)
    """
    
    deg_sign = {}   # gene -> sign
    detf_sign = {}  # tf -> sign (relies on tf_signs within loop)
    tf_avg = []  # collects the avg reg weights of TFs that regulates a DEG 
    deg_avg = [] # collects the reg weights of DEGs
    # as calculated by averaging all its DEG target's log2fc   
    
    for tf, edges in adjlist.items(): # avgs might not work

        tf_signs = [] # collects all tf's pos or neg regulation per deg
        tf_weight = []
        for edge in edges: # loops through each gene
            
            if len(edge) > 2: # if differential gene, then length > 2
                weight = edge[6] # log2fc col
                gene = edge[0] # ensembl id col
                sign = ''
                if weight > 0.0:
                    sign = "positive"
                else:
                    sign = "negative"
                
                deg_sign[gene] = sign # key = DEG, value = sign
                tf_signs.append(sign) # adds sign to list of tf reg
                deg_avg.append(weight) 
    
        tf_avg += sum(deg_avg) / len(deg_avg) # averages all deg weights per tf
        
        if tf_signs: # tf differential values are obtained by comparing counts 
  
            pos = tf_signs.count("positive")
            neg = tf_signs.count("negative")

            if pos > neg:
                detf_sign[tf] = "positive"
            elif neg > pos:
                detf_sign[tf] = "negative"
            else:
                detf_sign[tf] = "zero"

    summary = { 
        "DEGs_total": len(deg_sign),
        "DEGs_average":  sum(deg_avg) / len(deg_avg), # averages all degs
        "DEGs_positive": sum(1 for s in deg_sign.values() if s == "positive"),
        "DEGs_negative": sum(1 for s in deg_sign.values() if s == "negative"),
        
        # Below are regular TFS THAT REGULATE A DEG (so maybe special tfs)
        "TFs_total": len(detf_sign),
        "TFs_avg": sum(tf_avg) / len(tf_avg), # average TF
        "TFs_positive": sum(1 for s in detf_sign.values() if s == "positive"),
        "TFs_negative": sum(1 for s in detf_sign.values() if s == "negative"),
        }

    print(summary) # remove later
    return summary

def kmeans(adjlist): # struggles with clusters of different densities (aka evens out size)
    
    
    
    return

def calculate_clustering_coeff(adjlist): # input = adjlist like the one from A
    
    """
    Input:

    
    Output:  dict of (node:str, clustering_coeff:float)
    
    """
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