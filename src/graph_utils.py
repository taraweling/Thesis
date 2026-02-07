import random
import csv
import requests
import re # use regex
from collections import defaultdict
import numpy as np
# GRN in the form of adjacency list (adjlist) = dictionary:str (TFs) of lists:str (Genes) of tuples:(str,float) (Weights)
## eg {AHR:[(A1BG,0.47),(A1CF,-0.89)]..., AIRE:[(x,y),(a,b)]}

ENSEMBL_REST = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json"}

ENSEMBL_ID_RE = re.compile(r"^ENSG\d+$")
LOC_RE = re.compile(r"(chr)?(\w+):(\d+)-(\d+)")
XLOC_RE = re.compile(r"XLOC_\d+")

def gene2ensembl(query:str | None): # used chatgpt to complete loc_match portion of code 
    """
    Inputs a string containing a common name for a gene OR its location
    Outputs the ensembl ID if there, otherwise does not return
    """
    
    query = query.strip()

    # --- Case 1: Chromosomal location ---
    loc_match = re.match(r"(chr)?(\w+):(\d+)-(\d+)", query) 
    # Regex breakdown:
    # (chr)?        -> optional literal "chr" prefix
    # (\w+)         -> chromosome identifier (e.g. 1, 2, X, Y, MT)
    # :             -> literal colon
    # (\d+)         -> start coordinate
    # -             -> literal dash
    # (\d+)         -> end coordinate
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
    if re.match(r"XLOC_\d+", query):
        url = f"{ENSEMBL_REST}/xrefs/symbol/homo_sapiens/{query}"

        r = requests.get(url, headers=HEADERS)
        if not r.ok:
            return None

        xrefs = r.json()
        for x in xrefs:
            if x.get("type") == "gene":
                return x.get("id")
        return None
    # returned gene list is ordered by Ensembl’s internal logic (primarily coordinate order), not likelihood or biological relevance.
    
    # --- Case 3: Gene symbol / name ---
    url = f"{ENSEMBL_REST}/lookup/symbol/homo_sapiens/{query}"

    r = requests.get(url, headers=HEADERS)
    if not r.ok:
        return None

    data = r.json()
    return data.get("id")

def ensemblifylist(disorderlist):
    """
    Inputs:
        records: list of lists in format
            [GENEID, DISORDER, STUDY, YEAR, TISSUE, LOG2FC, PVAL]
    Outputs:
        same structure, but each GENEID converted to Ensembl ID if possible.
        If conversion fails, original GENEID is retained.
    """
    ensemblified = []

    for rec in disorderlist:
        gene = rec[0]
        ensembl_id = gene2ensembl(gene)
        ensembl_gene = ensembl_id if ensembl_id is not None else gene
        new_rec = [ensembl_gene] + rec[1:]
        ensemblified.append(new_rec)

    return ensemblified    

# underscores at start of fn name = internal helpers
def _extract_gene_strings(genedict):
    """
    Collect all *gene-like* strings appearing as dict keys or as first
    elements of tuple/list values. Non-gene metadata is ignored.
    """
    genes = set()

    for k, v in genedict.items():
        if isinstance(k, str):
            genes.add(k)

        # GRN-style: key -> list[(gene, weight)]
        if isinstance(v, list):
            for item in v:
                if isinstance(item, (tuple, list)) and isinstance(item[0], str):
                    genes.add(item[0])

        # DEG disorder_dict input rather than GRN: key is , value is metadata tuple = nothing to add
    return genes

def _bulk_symbol_lookup(symbols):
    """
    Resolve gene symbols → ENSG IDs via a single bulk POST.
    """
    if not symbols:
        return {}

    url = f"{ENSEMBL_REST}/lookup/symbol/homo_sapiens"
    r = requests.post(url, headers=HEADERS, json={"symbols": list(symbols)})
    if not r.ok:
        return {}

    data = r.json()
    return {
        sym: rec["id"]
        for sym, rec in data.items()
        if rec and "id" in rec
    }

def _resolve_special_cases(queries):
    """
    Handle location-based queries and XLOC identifiers individually.
    These cannot be bulk-resolved via the symbol endpoint.
    """
    resolved = {}

    for q in queries:
        # already ENSG
        if ENSEMBL_ID_RE.match(q):
            resolved[q] = q
            continue

        # genomic location
        m = LOC_RE.match(q)
        if m:
            chrom, start, end = m.group(2), m.group(3), m.group(4)
            url = f"{ENSEMBL_REST}/overlap/region/human/{chrom}:{start}-{end}"
            r = requests.get(url, headers=HEADERS, params={"feature": "gene"})
            if r.ok:
                genes = r.json()
                if genes:
                    resolved[q] = genes[0].get("id")
            continue

        # XLOC
        if XLOC_RE.match(q):
            url = f"{ENSEMBL_REST}/xrefs/symbol/homo_sapiens/{q}"
            r = requests.get(url, headers=HEADERS)
            if r.ok:
                for x in r.json():
                    if x.get("type") == "gene":
                        resolved[q] = x.get("id")
                        break

    return resolved

def ensemblify(genedict):
    """
    Efficiently convert all gene identifiers in a dictionary to ENSEMBL IDs.
    Uses one bulk POST for gene symbols, minimal individual calls for locations / XLOCs, zero calls for already-ENSEMBL IDs
    """

    # 1. collect all candidate gene strings
    genes = _extract_gene_strings(genedict)

    # 2. partition
    ensg = {g for g in genes if ENSEMBL_ID_RE.match(g)}
    locs = {g for g in genes if LOC_RE.match(g)}
    xlocs = {g for g in genes if XLOC_RE.match(g)}
    symbols = genes - ensg - locs - xlocs

    # 3. resolve
    mapping = {}
    mapping.update({g: g for g in ensg})
    mapping.update(_bulk_symbol_lookup(symbols))
    mapping.update(_resolve_special_cases(locs | xlocs))

    # 4. rebuild dictionary with ENSG keys/values
    newdict = {}

    for k, v in genedict.items(): # changed from newdict to genedict
        new_k = mapping.get(k, k)

        if isinstance(v, list):
            new_v = []
            for gene, weight in v:
                new_gene = mapping.get(gene, gene)
                new_v.append((new_gene, weight))
            newdict[new_k] = new_v
        else:
            newdict[new_k] = v

    return newdict

def oldensemblify(deggrn): # inputs dict containing GRN to entirely ensemblify all genes. SPEED UP VIA BULK CALL - am i making unecessary calls? 
    for k in list(deggrn): # loops through each tf and its values
        new = gu.gene2ensembl(k)
        deggrn[new] = deggrn.pop(k) # ensemblfies all keys
        v = deggrn[new] # looks at each TF's targets
        for i, gene in enumerate(v): #
            new = gu.gene2ensembl(gene[0])
            newgene = new if new is not None else gene[0] # if input is already ensembl, replace 'None' output with original input
            v[i] = (newgene,gene[1])
    return deggrn

def adjlist2adjmat(adjlist:dict): # works for weighted adjlist. if including adjmat as a file, include output location
    
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

def adjlist2edgelist(adjlist:dict):
    
    """
    Convert a directed adjacency list into a directed edge list,
    preserving all edge-associated metadata.

    Accepts BOTH formats:

    (1) Simple GRN:
        {
            TF: [(Gene, weight), ...]
        }

    (2) DEG-annotated GRN:
        {
            TF: [
                (Gene, weight, disorder, study, year, tissue, log2fc, pval),
                ...
            ]
        }

    Output:
        [
            [TF, Gene, weight],
            ...
        ]
        OR
        [
            [TF, Gene, weight, disorder, study, year, tissue, log2fc, pval],
            ...
        ]
    """

    edgelist = []

    for tf, edges in adjlist.items():
        for edge in edges:

            # minimal case: (Gene, weight)
            if len(edge) == 2:
                gene, weight = edge
                edgelist.append([tf, gene, weight])

            # annotated case:
            else:
                gene = edge[0]
                metadata = edge[1:]
                edgelist.append([tf, gene, *metadata])

    
    return edgelist

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

def disorder_list(filename:str, *disorder:str): # generates the joined list of list of n number of disorders, fixed using chatgpt
    """
    Inputs specific disorder(s) names and the name of csv file storing all disorders' DEGs in the following format:
    DISORDER | STUDY | YEAR | TISSUE | GENEID | LOG2FC | PVAL
    
    Outputs a list of lists of that disorder/set of disorder's DEGs with the following format:
    [[GENE1, DISORDER1, STUDY1, YEAR1, TISSUE1, LOG2FC1, PVAL1],[GENE2, DISORDER2, STUDY2, YEAR2, TISSUE2, LOG2FC2, PVAL2],[GENE3....]]
    
    Output preserves all duplicates for the same gene.
    """
    disorder_set = set(d.strip() for d in disorder)
    records = []

    with open(filename, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            dname = row["DISORDER"].strip()
            if dname not in disorder_set:
                continue

            record = [
                row["GENEID"].strip(),
                dname,
                row["STUDY"].strip(),
                int(row["YEAR"]),
                row["TISSUE"].strip(),
                float(row["LOG2FC"]),
                float(row["PVAL"]),
            ]
            records.append(record)

    return records

def deg_grn_tfsonly(grn,disorderlist): 
    """
    Merge information from the GRN adjacency list (dict) and disorders list of lists such that the resulting dictionary looks as follows: 
    {TF1:[(Gene1, GRNedgeweight1, disorder1, study1, year1, tissue1, log2fc1, pval1), 
    (Gene2, GRNedgeweight2, disorder2, study2, year2, tissue2, log2fc2, pval2)], TF2:[(etc)]}
    
    In the example above, the keys can hold lists with repeat tuples, 
    eg TF2:[(Gene1, GRNedgeweight1, disorder1, study1, year1, tissue1, log2fc1, pval1),
    (Gene3, GRNedgeweight3, disorder3, study3, year3, tissue3, log2fc3, pval3)]
    
    """
    gene_to_disorders = {}
    for rec in disorderlist:
        gene_id = rec[0]
        info = tuple(rec[1:])  # (DISORDER, STUDY, YEAR, TISSUE, LOG2FC, PVAL)
        gene_to_disorders.setdefault(gene_id, []).append(info)

    disorder_genes = set(gene_to_disorders.keys())

    finalgrn = {}

    for tf, targets in grn.items():

        # ONLY restrict by TF presence in disorderlist
        if tf not in disorder_genes:
            continue

        merged = []

        for gene, weight in targets:

            # If gene has disorder info, attach all of it
            if gene in gene_to_disorders:
                for disorder_info in gene_to_disorders[gene]:
                    merged.append((gene, weight) + disorder_info)

            # Otherwise keep gene + weight only
            else:
                merged.append((gene, weight))

        finalgrn[tf] = merged
    return finalgrn

def deg_grn_both(grn, disorderlist):

    gene_to_disorders = {}

    for rec in disorderlist:
        gene_to_disorders.setdefault(rec[0], []).append(tuple(rec[1:]))

    degset = set(gene_to_disorders)

    finalgrn = {}

    for tf, targets in grn.items():

        if tf not in degset:
            continue

        merged = []

        for edge in targets:
            gene = edge[0]
            weight = edge[1]

            if gene not in gene_to_disorders:
                continue

            for disorder_info in gene_to_disorders[gene]:
                merged.append((gene, weight) + disorder_info)

        if merged:
            finalgrn[tf] = merged

    return finalgrn


def expand_node_attributes(adjlist:dict, genes:list): # UNFINISHED
    """
    
    
    """
    return