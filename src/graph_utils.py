import random
import csv
import requests
import math
import re # use regex
from collections import defaultdict
import numpy as np
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
# GRN in the form of adjacency list (adjlist) = dictionary:str (TFs) of lists:str (Genes) of tuples:(str,float) (Weights)
## eg {AHR:[(A1BG,0.47),(A1CF, 0.89)]..., AIRE:[(x,y),(a,b)]}

ENSEMBL_REST = "https://rest.ensembl.org"
HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",}

# Shared session with retries to reduce transient Ensembl timeouts.
_SESSION = None
_DEFAULT_TIMEOUT = (3.05, 10)

_RETRY = Retry(
    total=3,
    connect=3,
    read=3,
    backoff_factor=0.5,
    status_forcelist=(429, 500, 502, 503, 504),
    allowed_methods=("GET", "POST"),
    raise_on_status=False,
)


def _get_session():
    global _SESSION
    if _SESSION is None:
        s = requests.Session()
        adapter = HTTPAdapter(max_retries=_RETRY)
        s.mount("https://", adapter)
        s.mount("http://", adapter)
        _SESSION = s
    return _SESSION


def _request(method, url, **kwargs):
    timeout = kwargs.pop("timeout", _DEFAULT_TIMEOUT)
    try:
        return _get_session().request(method, url, timeout=timeout, **kwargs)
    except requests.exceptions.RequestException:
        return None

ENSEMBL_ID_RE = re.compile(r"^ENSG\d+$")
LOC_RE = re.compile(r"(chr)?(\w+):(\d+) (\d+)") # regex
XLOC_RE = re.compile(r"XLOC_\d+")

# purpose of gene2ensembl vs ensemblify: single lookups, debugging, 
def gene2ensembl(query:str | None): # used chatgpt to complete loc_match portion of code 
    """
    Inputs a string containing a common gene name OR a genomic location.
    Returns the Ensembl gene ID if found; otherwise returns None.
    """
    
    query = query.strip()

    #     Case 1: Chromosomal location    
    loc_match = re.match(r"(chr)?(\w+):(\d+) (\d+)", query) 
    # Regex breakdown:
    # (chr)?         > optional literal "chr" prefix
    # (\w+)          > chromosome identifier (e.g. 1, 2, X, Y, MT)
    # :              > literal colon
    # (\d+)          > start coordinate
    #                > literal dash
    # (\d+)          > end coordinate
    if loc_match:
        chrom = loc_match.group(2)
        start = loc_match.group(3)
        end = loc_match.group(4)

        url = f"{ENSEMBL_REST}/overlap/region/human/{chrom}:{start} {end}"
        params = {"feature": "gene"}

        r = _request("GET", url, headers=HEADERS, params=params, timeout=10)
        if r is None or not r.ok:
            return None

        genes = r.json()
        if not genes:
            return None

        # Return the first overlapping gene  > does this mean the most likely match?
        return genes[0].get("id")
    
    #     Case 2: XLOC    
    if re.match(r"XLOC_\d+", query):
        url = f"{ENSEMBL_REST}/xrefs/symbol/homo_sapiens/{query}"

        r = _request("GET", url, headers=HEADERS, timeout=10)
        if r is None or not r.ok:
            return None

        xrefs = r.json()
        for x in xrefs:
            if x.get("type") == "gene":
                return x.get("id")
        return None
    # returned gene list is ordered by Ensembl’s internal logic (primarily coordinate order), not likelihood or biological relevance.
    
    #     Case 3: Gene symbol / name    
    url = f"{ENSEMBL_REST}/lookup/symbol/homo_sapiens/{query}"

    r = _request("GET", url, headers=HEADERS, timeout=10)
    if r is None or not r.ok:
        return None

    data = r.json()
    return data.get("id")
        
def ensembl2gene(query:str | None):
    # returns a gene symbol (hgnc) for any ensembl sequence. use for visualizing
    # uses the shared the HGNC label cache across the 
    # per‑disorder visualizations so we make far fewer requests.
    
    if query is None:
        return None

    query = query.strip()
    if not query:
        return None

    # Case 1: Chromosomal location
    loc_match = re.match(r"(chr)?(\w+):(\d+) (\d+)", query)
    if loc_match:
        chrom = loc_match.group(2)
        start = loc_match.group(3)
        end = loc_match.group(4)

        url = f"{ENSEMBL_REST}/overlap/region/human/{chrom}:{start} {end}"
        params = {"feature": "gene"}

        r = _request("GET", url, headers=HEADERS, params=params, timeout=10)
        if r is None or not r.ok:
            return None

        genes = r.json()
        if not genes:
            return None

        # Return the first overlapping gene name if available
        return genes[0].get("external_name") or genes[0].get("id")

    # Case 2: XLOC
    if re.match(r"XLOC_\d+", query):
        url = f"{ENSEMBL_REST}/xrefs/symbol/homo_sapiens/{query}"

        r = _request("GET", url, headers=HEADERS, timeout=10)
        if r is None or not r.ok:
            return None

        xrefs = r.json()
        for x in xrefs:
            if x.get("type") == "gene":
                gene_id = x.get("id")
                if not gene_id:
                    continue
                lookup = _request("GET", f"{ENSEMBL_REST}/lookup/id/{gene_id}", headers=HEADERS, timeout=10)
                if lookup is None or not lookup.ok:
                    return gene_id
                data = lookup.json()
                return data.get("display_name") or data.get("external_name") or gene_id
        return None

    # Case 3: Ensembl ID (or other stable ID)
    url = f"{ENSEMBL_REST}/lookup/id/{query}"

    r = _request("GET", url, headers=HEADERS, timeout=10)
    if r is None or not r.ok:
        return None

    data = r.json()
    return data.get("display_name") or data.get("external_name") or data.get("id")

# underscores at start of fn name = internal helpers

def _extract_gene_strings(genedict):
    """
    Collect all *gene-like* strings found as:
    1) dictionary keys, and
    2) the first element of tuple/list values.
    Any non-gene metadata is ignored.
    """
    genes = set()

    for k, v in genedict.items():
        if isinstance(k, str):
            genes.add(k)

        # GRN style: key  > list[(gene, weight)]
        if isinstance(v, list):
            for item in v:
                if isinstance(item, (tuple, list)) and isinstance(item[0], str):
                    genes.add(item[0])

        # DEG disorder_dict input rather than GRN: key is , value is metadata tuple = nothing to add
    return genes

def overlap_summary_row(name, degset, grn_tfs, grn_nodes):
    """
    Build a per-disorder overlap row for easy CSV comparison.
    Inputs are sets (or iterables) of DEG genes, GRN TFs, and all GRN nodes.
    """
    degset = set(degset) if degset is not None else set()
    grn_tfs = set(grn_tfs) if grn_tfs is not None else set()
    grn_nodes = set(grn_nodes) if grn_nodes is not None else set()

    deg_count = len(degset)
    grn_tf_count = len(grn_tfs)
    grn_node_count = len(grn_nodes)

    overlap_tfs = len(grn_tfs & degset)
    overlap_nodes = len(grn_nodes & degset)

    def _ratio(num, den):
        return round(num / den, 6) if den else 0

    return {
        "disorder": name,
        "deg_count": deg_count,
        "grn_tf_count": grn_tf_count,
        "grn_node_count": grn_node_count,
        "deg_overlap_tfs": overlap_tfs,
        "deg_overlap_grn_nodes": overlap_nodes,
        "pct_degs_that_are_tfs": _ratio(overlap_tfs, deg_count),
        "pct_degs_in_grn_nodes": _ratio(overlap_nodes, deg_count),
        "pct_grn_tfs_in_degs": _ratio(overlap_tfs, grn_tf_count),
        "pct_grn_nodes_in_degs": _ratio(overlap_nodes, grn_node_count),
    }

def write_overlap_summary(rows, outfile):
    
    #Write overlap summary rows to a CSV with stable, comparable columns.
    
    if not rows:
        return

    fieldnames = [
        "disorder",
        "deg_count",
        "grn_tf_count",
        "grn_node_count",
        "deg_overlap_tfs",
        "deg_overlap_grn_nodes",
        "pct_degs_that_are_tfs",
        "pct_degs_in_grn_nodes",
        "pct_grn_tfs_in_degs",
        "pct_grn_nodes_in_degs",
    ]

    with open(outfile, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

def _bulk_symbol_lookup(symbols):
    # Resolve gene symbols → ENSG IDs using one bulk POST call.
    if not symbols:
        return {}

    url = f"{ENSEMBL_REST}/lookup/symbol/homo_sapiens"
    r = _request("POST", url, headers=HEADERS, json={"symbols": list(symbols)}, timeout=20)
    if r is None or not r.ok:
        return {}

    data = r.json()
    return {
        sym: rec["id"]
        for sym, rec in data.items()
        if rec and "id" in rec}

def _resolve_special_cases(queries):
    
    #Handle location-based queries and XLOC identifiers individually.
    #These cannot be bulk-resolved via the Ensembl batch endpoint.
    
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
            url = f"{ENSEMBL_REST}/overlap/region/human/{chrom}:{start} {end}"
            r = _request("GET", url, headers=HEADERS, params={"feature": "gene"}, timeout=10)
            if r is not None and r.ok:
                genes = r.json()
                if genes:
                    resolved[q] = genes[0].get("id")
            continue

        # XLOC
        if XLOC_RE.match(q):
            url = f"{ENSEMBL_REST}/xrefs/symbol/homo_sapiens/{q}"
            r = _request("GET", url, headers=HEADERS, timeout=10)
            if r is not None and r.ok:
                for x in r.json():
                    if x.get("type") == "gene":
                        resolved[q] = x.get("id")
                        break

    return resolved

# purpose of ensemblify and ensemblifylist: thousands of genes, DEGs, GRNs
def ensemblify(genedict):
    """
    Convert all gene names/symbols in a dictionary to ENSEMBL IDs.
    1) Use a single bulk POST for gene symbols.
    2) Use minimal individual calls for locations/XLOCs.
    3) Make zero calls for already-ENSEMBL identifiers.
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

def ensemblifylist(disorderlist):
    """
    Inputs: records: list of lists in format [GENEID, DISORDER, STUDY, YEAR, TISSUE, LOG2FC, PVAL]
    
    Outputs: same structure, but each GENEID converted to Ensembl ID if possible.
    If conversion fails, original GENEID is retained.

    Conversion is done in bulk using helper functions.
    """
    genes = {rec[0] for rec in disorderlist} # collect unique genes
    
    # partition by type of gene ID
    ensg = {g for g in genes if ENSEMBL_ID_RE.match(g)}
    locs = {g for g in genes if LOC_RE.match(g)}
    xlocs = {g for g in genes if XLOC_RE.match(g)}
    symbols = genes - ensg - locs - xlocs
    
    # make names consistent
    mapping = {}

    # already ENSG
    mapping.update({g: g for g in ensg})

    # bulk symbol  > ENSG
    if symbols:
        mapping.update(_bulk_symbol_lookup(symbols))

    # individual chrom coords or XLOC
    if locs or xlocs:
        mapping.update(_resolve_special_cases(locs | xlocs))
        
    """
    Previous per-record implementation (kept for reference):
    ensemblified = []

    for rec in disorderlist:
        gene = rec[0]
        ensembl_id = gene2ensembl(gene)
        ensembl_gene = ensembl_id if ensembl_id is not None else gene
        new_rec = [ensembl_gene] + rec[1:]
        ensemblified.append(new_rec)

    return ensemblified
    """
    
    # rebuild list of DEGs
    out = []

    for rec in disorderlist:
        gene = rec[0] # extracts 1st column = gene identifier
        out.append([mapping.get(gene, gene), *rec[1:]]) # references gene in mapping dict
        #if found, returns new val otherwise returns original gene. also skips first column

    return out

def oldensemblify(deggrn): # inputs dict containing GRN to entirely ensemblify all genes. SPEED UP VIA BULK CALL   am i making unecessary calls? 
    for k in list(deggrn): # loops through each tf and its values
        new = gene2ensembl(k)
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
    1) Simple GRN:
       { TF: [(Gene, weight), ...] }

    2) DEG annotated GRN:
       {
         TF: [
           (Gene, weight, disorder, study, year, tissue, log2fc, pval),
           ... ]}

    Output:
      [[TF, Gene, weight], ...]
      OR
      [[TF, Gene, weight, disorder, study, year, tissue, log2fc, pval], ...]
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
    Inputs: any number of adjacency lists containing GRNs.
    Behavior:
      - Merge repeated TF→gene pairs across graphs.
      - Average duplicate weights using numpy.
    Output: an averaged GRN (biological validity not established).

    Biologically, this is a naive consensus step across multiple inferred GRNs.
    Averaging edge weights is one simple way to aggregate evidence across
    networks or cohorts, but it does not preserve individual-specific
    regulatory variation highlighted in GRN-based disease studies. For
    example, network approaches (e.g., PANDA/LIONESS) emphasize that
    regulatory edges can differ across individuals or conditions and may
    reveal dysregulation not captured by expression alone. (cite: anwer_transcriptional_2025)
    """
    edge_weights = {}  # tracks repeats across grns for each tf gene pair, eg {(AHR,A1BG):[weight in GRN1, weight in GRN2...]}
    

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
    Inputs:
      - oldgraph: adjacency list for a complete GRN
      - disorder: list of TFs to keep
    Output:
      - filtered adjacency list with only those TF keys
    """
    newgraph = {}
    
    for tf, edges in oldgraph.items(): # edges = (gene,weight)
        if tf in disorder and tf not in newgraph:
            newgraph[tf] = edges            
    
    return newgraph
    
def make_adjlist(filename:str, threshold:float): 
    """
    Inputs:
      - filename: path to a CSV matrix where first row is target genes
      - threshold: minimum weight to keep an edge
    Output:
      - adjacency list {TF: [(Gene, Weight), ...]}
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

def disorder_list(filename:str, *disorder:str): # generates the concatenated list of list of n number of disorders, fixed using chatgpt
    """
    Inputs:
      - disorder names
      - CSV with columns:
        DISORDER | STUDY | YEAR | TISSUE | GENEID | LOG2FC | PVAL
    Output:
      - list of records:
        [GENE, DISORDER, STUDY, YEAR, TISSUE, LOG2FC, PVAL]
    Duplicates are preserved.
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
                float(row["PVAL"]),]
            records.append(record)

    return records

def de_grn_tfsonly(grn,disorderlist): # aka enrich GRN with DEG info
    """
    Merge GRN adjacency list with disorder metadata, producing:
    {TF: [(Gene, GRN_weight, disorder, study, year, tissue, log2fc, pval), ...]}

    This version keeps only TFs that are themselves DEGs and only targets
    that are also DEGs. TF keys can still contain repeated tuples.

    Biological rationale:
    - TF networks are central to maintaining transcriptional states, and only
      a subset of TF perturbations may drive major transcriptomic shifts.
      (cite: nishiyama_systematic_2013)
    - By intersecting DE TFs with DE targets, this focuses on edges most
      plausibly connected to disease-associated expression changes.
    - Direction of regulation (activation vs repression) can be context-
      dependent, so enrichment alone does not imply a single regulatory sign.
      (cite: martinez-corral_emergence_2024)
    """
    gene_to_disorders = {}

    """
    Build a lookup from gene → disorder metadata so GRN edges can be
    annotated with expression context; this ties network structure to
    disease-associated transcriptomic shifts. (cite: anwer_transcriptional_2025)
    """
    for rec in disorderlist:
        gene_id = rec[0]
        info = tuple(rec[1:])
        gene_to_disorders.setdefault(gene_id, []).append(info)

    degset = set(gene_to_disorders.keys())

    finalgrn = {}

    for tf, targets in grn.items():

        """
        Keep only TFs that are themselves DEGs, reflecting the idea that
        a subset of TFs drive the strongest transcriptomic changes.
        (cite: nishiyama_systematic_2013)
        """
        if tf not in degset:
            continue

        merged = []

        for gene, weight in targets:

            """
            Only keep targets that are DEGs; attach all available disorder
            rows because TF effects can vary by context.
            (cite: martinez-corral_emergence_2024)
            """
            if gene not in degset:
                continue

            for disorder_info in gene_to_disorders[gene]:
                merged.append((gene, weight) + disorder_info)

        if merged:
            finalgrn[tf] = merged

    return finalgrn

def de_grn_both(grn, disorderlist):
    
    """
    Merge GRN adjacency list with disorder metadata, keeping only
    TFs and targets that are DEGs. Attaches ALL disorder rows per gene.

    Biological rationale:
    - TF-centric GRN analysis highlights a small subset of regulators with
      outsized influence on expression programs. (cite: nishiyama_systematic_2013)
    - Keeping all disorder rows preserves multiple observations of the same
      gene, which is important because TF effects can vary by context and
      may show mixed activation/repression across conditions. (cite: martinez-corral_emergence_2024)
    """

    # Build gene from list of disorder metadata
    gene_to_disorders = {}

    """
    Preserve all observations of a gene across studies/tissues because
    regulatory effects can depend on molecular context. (cite: martinez-corral_emergence_2024)
    """
    for rec in disorderlist:
        gene_id = rec[0]
        info = tuple(rec[1:])  # (disorder, study, year, tissue, log2fc, pval)
        gene_to_disorders.setdefault(gene_id, []).append(info)

    degset = set(gene_to_disorders.keys())

    finalgrn = {}

    for tf, targets in grn.items():

        """
        Keep only TFs in the DEG set to focus on regulators directly
        associated with expression shifts. (cite: nishiyama_systematic_2013)
        """
        if tf not in degset:
            continue

        merged = []

        for gene, weight in targets:

            if gene not in degset:
                continue

            # attach ALL disorder rows for this gene
            for disorder_info in gene_to_disorders[gene]:
                merged.append((gene, weight) + disorder_info)

        if merged:
            finalgrn[tf] = merged

    return finalgrn

def de_grn_degsonly(grn, disorderlist):

    """
    Identify regulatory edges where:
      - TF is NOT a DEG
      - target gene IS a DEG
    Output format:
      {TF:[(Gene, GRNedgeweight, disorder, study, year, tissue, log2fc, pval), ...]}

    Biological rationale:
    - Differential expression in targets can be driven by upstream regulators
      even if the TF itself is not differentially expressed, consistent with
      GRN-based disease models where regulatory edges shift without large
      TF expression changes. (cite: anwer_transcriptional_2025)
    - This separation helps isolate potential regulatory influence that is
      not captured by TF differential expression alone. (cite: nishiyama_systematic_2013)
    """

    # Build gene  to disorder "metadata" dictionary
    gene_to_disorders = {}

    """
    Map genes to all disorder observations so downstream edges can reflect
    multi-study context rather than a single measurement. (cite: anwer_transcriptional_2025)
    """
    for rec in disorderlist:
        gene_id = rec[0]
        info = tuple(rec[1:])
        gene_to_disorders.setdefault(gene_id, []).append(info)

    degset = set(gene_to_disorders.keys())

    finalgrn = {}

    for tf, targets in grn.items():

        # Skip TFs that ARE DEGs
        if tf in degset:
            continue

        merged = []

        for gene, weight in targets:

            # Only keep targets that ARE DEGs
            if gene not in degset:
                continue

            # Attach all disorder metadata rows
            for disorder_info in gene_to_disorders[gene]:
                merged.append((gene, weight) + disorder_info)

        if merged:
            finalgrn[tf] = merged

    return finalgrn

def _tf_targets(adj):
    """
    Build TF -> set(target genes) from an adjacency list.
    Supports simple (gene, weight) edges and annotated edges where the gene
    identifier remains in index 0.
    """
    tf_targets = {}

    for tf, edges in adj.items():
        tf_targets[tf] = {
            edge[0]
            for edge in edges
            if isinstance(edge, (list, tuple)) and len(edge) > 0
        }

    return tf_targets

def feed_forward_loops(adj, degset):
    """
    FFL definition:
      TF A regulates TF B
      TF A regulates DEG C
      TF B regulates DEG C
    """
    loops = []

    tf_targets = _tf_targets(adj)
    tfs = set(adj.keys())
    degset = set(degset) if degset is not None else set()

    for tfA in tfs:
        targetsA = tf_targets.get(tfA, set())

        # TF B must also be a TF node.
        tf_targets_in_A = targetsA & tfs
        for tfB in tf_targets_in_A:
            targetsB = tf_targets.get(tfB, set())
            shared_targets = targetsA & targetsB & degset

            for gene in shared_targets:
                loops.append({
                    "TF_A": tfA,
                    "TF_B": tfB,
                    "target_DEG": gene,
                })

    return loops

def feedback_loops(adj):
    """
    Feedback loop definition:
      TF A regulates TF B
      TF B regulates TF A
    """
    loops = []

    tf_targets = _tf_targets(adj)
    tfs = set(adj.keys())

    for tfA in tfs:
        targetsA = tf_targets.get(tfA, set())
        tf_targets_in_A = targetsA & tfs

        for tfB in tf_targets_in_A:
            targetsB = tf_targets.get(tfB, set())
            if tfA in targetsB:
                loops.append({
                    "TF_A": tfA,
                    "TF_B": tfB,
                })

    return loops

def detect_regulatory_loops(degs, grn_by_name, graph_label="grn"):
    """
    Run per-disorder feed-forward and feedback loop detection over a GRN map.
    Returns:
      (ffl_by_name, fbl_by_name, ffl_rows, fbl_rows)
    """
    ffl_by_name = {}
    fbl_by_name = {}
    ffl_rows = []
    fbl_rows = []

    print(f"\nsearching for regulatory loops ({graph_label})")

    for name, data in degs.items():
        if name not in grn_by_name:
            continue
        if not data:
            continue

        print("\nloop detection for", name, f"({graph_label})")

        adj = grn_by_name[name]
        degset = {row[0] for row in data if row}

        ffls = feed_forward_loops(adj, degset)
        fbls = feedback_loops(adj)

        ffl_by_name[name] = ffls
        fbl_by_name[name] = fbls

        for rec in ffls:
            ffl_rows.append({
                "graph": graph_label,
                "disorder": name,
                "TF_A": rec["TF_A"],
                "TF_B": rec["TF_B"],
                "target_DEG": rec["target_DEG"],
            })

        for rec in fbls:
            fbl_rows.append({
                "graph": graph_label,
                "disorder": name,
                "TF_A": rec["TF_A"],
                "TF_B": rec["TF_B"],
            })

        print(
            name,
            "FFLs:",
            len(ffls),
            "FBLs:",
            len(fbls),
        )

    return ffl_by_name, fbl_by_name, ffl_rows, fbl_rows

def expand_node_attributes(adjlist:dict, genes:list): # UNFINISHED
    """
    Placeholder for future node attribute expansion logic.
    """
    return
