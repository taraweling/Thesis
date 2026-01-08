import random
import csv
import requests
import re # use regex
from collections import defaultdict # cool 
# GRN =  dictionary:str (TFs) of lists:str (Genes) of tuples:(str,float) (Weights)
## eg {AHR:[(A1BG,0.47),(A1CF,-0.89)]..., AIRE:[(x,y),(a,b)]}
ENSEMBL_REST = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json"}


def gene2ensembl(query:str | None): # used chatgpt to generate code for case 1 + debug code
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

def merge_adjlist(*oldgraph): # should be able to input any num of adjlists to merge
    """
    Inputs n number of adjacency list containing a GRN, combine all
    Merge repeated TF-gene pairs between graphs, averaging the weights for duplicates
    
    Outputs an averaged GRN, with the caveat that this may have limited biological applicability as it is previously untested
    """
    newgraph = {}
    seen_edges = set() # tracks repeats
    
    for graph in oldgraph: # loops through each GRN 
        for tf, edges in graph.items(): # loops through both key (tf) and value (tuple of (gene,weight))
            if tf not in newgraph:
                newgraph[tf] = []
            
            for gene, weight in edges:
                edge = (tf, gene)
                print(weight)
                for i in newgraph:
                    print(type(i))
                    if edge in i:
                        print(gene)
                if edge not in seen_edges: 
                    newgraph[tf].append((gene,weight))
                    seen_edges.add(edge)  
    
    return newgraph

def expand_node_attributes(adjlist:dict, genes:list):
    """
    
    
    """
    return

def degs_disorder(filename:str, *disorder:str):
    """
    Inputs specific disorder(s) and name of file with all disorders's DEGs in the following format:
    DISORDER | STUDY | YEAR | TISSUE | GENEID | LOG2FC | PVAL
    
    Outputs a list of that disorder/set of disorder's DEGs
    """
    degs = []
    
    return degs
     
def filter_adjlist(adjlist:dict, disorder:list): # should be *disorder to take mult inputs?
    """
    Inputs an adjacency list representing a complete GRN
    Outputs filtered adjlist reducing the keys in a GRN to just those of the inputted list
    """
    
    
    return 
    

def make_adjlist(filename:str, threshold:float): 
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

