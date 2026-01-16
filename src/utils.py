import random
import csv
import requests
import re # use regex
import numpy as np
# GRN in the form of adjacency list (adjlist) = dictionary:str (TFs) of lists:str (Genes) of tuples:(str,float) (Weights)
## eg {AHR:[(A1BG,0.47),(A1CF,-0.89)]..., AIRE:[(x,y),(a,b)]}
ENSEMBL_REST = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json"}


def gene2ensembl(query:str | None): # used chatgpt to generate code 
    """
    Inputs a string containing a common name for a gene OR its location
    Outputs the ensembl ID if there, otherwise does not return
    """
    
    query = query.strip()

    # --- Case 1: Chromosomal location ---
    loc_match = re.match(r"(chr)?(\w+):(\d+)-(\d+)", query) # what does this do
    if loc_match:
        chrom = loc_match.group(2)
        start = loc_match.group(3)
        end = loc_match.group(4)

        url = f"{ENSEMBL_REST}/overlap/region/human/{chrom}:{start}-{end}"
        params = {"feature": "gene"}

        r = requests.get(url, headers=HEADERS, params=params)
        if not r.ok:
            return None

        genes = r.json()
        if not genes:
            return None

        # Return the first overlapping gene -> does this mean the most likely match?
        return genes[0].get("id")
    # --- Case 2: XLOC ---
    
    
    # --- Case 3: Gene symbol / name ---
    url = f"{ENSEMBL_REST}/lookup/symbol/homo_sapiens/{query}"

    r = requests.get(url, headers=HEADERS)
    if not r.ok:
        return None

    data = r.json()
    return data.get("id")

def adjlist2adjmat(adjlist:dict): # if including adjmat as a file, include output location
    
    
    return

def merge_adjlist(*oldgraph):
    """
    Inputs n number of adjacency list containing a GRN, combine all
    Merge repeated TF-gene pairs between graphs, averaging the weights for duplicates using numpy
    
    Outputs an averaged GRN, with the caveat that this may have limited biological applicability as it is previously untested
    """
    edge_weights = {}  # tracks repeats across grns for each tf-gene pair, eg {(AHR,A1BG):[weight in GRN1, weight in GRN2...]}
    

    for graph in oldgraph:
        for tf, edges in graph.items():
            for gene, weight in edges:
                edge = (tf, gene)
                if edge not in edge_weights:
                    edge_weights[edge] = []
                edge_weights[edge].append(weight) # collects repeat weights

    newgraph = {}
    for (tf, gene), weights in edge_weights.items():
        mean_weight = float(np.mean(weights))
        if tf not in newgraph:
            newgraph[tf] = []
        newgraph[tf].append((gene, mean_weight))

    return newgraph

def old_merge_adjlist(*oldgraph): # DOES NOT AVERAGE REPEATS, CHOOSES FIRST WEIGHT
    
    newgraph = {}
    seen_edges = set() # tracks repeats, not necessary since dicts already remove dups
    
    for graph in oldgraph: # loops through each GRN 
        for tf, edges in graph.items(): # loops through both key (tf) and value (tuple of (gene,weight))
            if tf not in newgraph:
                newgraph[tf] = []
            
            for gene, weight in edges:
                edge = (tf, gene)

                if edge not in seen_edges: 
                    newgraph[tf].append((gene,weight))
                    seen_edges.add(edge)  
    
    return newgraph

def expand_node_attributes(adjlist:dict, genes:list): # UNFINISHED
    """
    
    
    """
    return

     
def filter_adjlist(oldgraph:dict, disorder:list): # should be *disorder to take mult inputs?
    """
    Inputs an adjacency list representing a complete GRN
    Outputs filtered adjlist reducing the keys in a GRN to just those of the inputted list
    """
    newgraph = {}
    
    for tf, edges in oldgraph.items(): # edges = (gene,weight)
        if tf in disorder and tf not in newgraph:
            newgraph[tf] = edges            
    
    return newgraph
    
def make_adjlist(filename:str, threshold:float): # made with help from chatgpt
    """
    Inputs path:str to csv of genes (str); data has formatting requirements and a weight for the necessary threshold
    
    Outputs a dictionary:str (TFs) of lists:str (Genes) of tuples:(str,float) (Weights)
    """
    processed_data = [] # list to hold processed data
    adj = {}
    
    with open(filename, "r", newline="") as file:
        reader = csv.reader(file)
        header = next(reader) # first row: gene names
        target_genes = header[1:] # skip TF column

        for row in reader:
            tf = row[0]
            values = row[1:]

            for gene, val in zip(target_genes, values): 
                try:
                    weight = float(val)
                except (ValueError, TypeError):
                    continue # skip missing or invalid values
                
                if weight > threshold:
                    adj.setdefault(tf, [])
                    adj[tf].append((gene, weight))
                    
    return adj #adjacency list of the graph

def viz_graph(adjlist): # UNFINISHED #saves an html file of visualized graph
    
    
    
    
    return


def kmeans(adjlist): # UNFINISHED # struggles with clusters of different densities (aka evens out size)
    
    
    
    return