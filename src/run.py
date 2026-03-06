import csv
from pyvis.network import Network 
import networkx as nx
import matplotlib.pyplot as plt
import graph_utils as gu
import graph_viz as gv 
import graph_algos as ga 
import json
import py4cytoscape as pyc
import pandas as pd
import time # start = time.time(); multi_time = time.time() - start; print(f"Multi-threaded:  {multi_time:.3f}s")

# INSTRUCTIONS: RUN WHILE IN SRC FOLDER
def main():
    
    # obtain adjlist for each disorder (handle in main or helper fn?)
    ## what information do I want in my graph?
    """
    as a bipartite graph, node attribute should either be TF or gene. 
    the TFs are DEGs and the target genes can be either DEGs or not-differentially expressed genes in the GRAND regulatory network 
    """
    
    # test make_adjlist
    #testadjlist = utils.make_adjlist('data/test.csv',0) # truong paper uses edge weight threshold of 0.7. 
    # would this change since I'm merging 2 adjlists
    #print(testadjlist)
    
    # get the combined GRAND GRN (what check can I use?)
    brainother = gu.make_adjlist('data/Brain_Other.csv',0.7) 
    brainbg = gu.make_adjlist('data/Brain_Basal_Ganglia.csv',0.7)

    # test merge_adjlist's weight-averaging feature: AHR-A1BG SHOULD have a weight of ~0.2
    #testadjlist2 = {'AHR': [('A1BG', 0.0), ('A4GALT', 0.314864), ('AHR', 0.131604), ('ALX3', 0.208409)]}
    #testmerge = utils.merge_adjlist(testadjlist,testadjlist2)
    #print(testmerge,'\n')

    # merge the two brain GRNs 
    brains= gu.merge_adjlist(brainother, brainbg) # NOT WORKING 
    grn = gu.ensemblify(brains)
    
    """Get the first key
    first_node = next(iter(brains))
    print("non ensemble",first_node)"""

    # should I do a bunch of stats calls on the merged brain grn to compare to individual?
    
    # get adjlists per disorder by inputting the name of the disorder
    ## (options = AD, ADHD, BD, SZ, MDD, OCD) and the file location, outputting a list of lists   
    test = gu.disorder_list('data/DEGDataSample.csv','BD')

    # convert geneIDs in disorders list of lists to ENSEMBL IDs
    degs = gu.ensemblifylist(test)
    

    # reduce GRN to just keys from disorders (REPLACED WITH DEG_GRN_BOTH VS DEG_GRN_TFS_ONLY)
    #deggrn = gu.filter_adjlist(brain,disorders)
    #print(len(grn),len(grn.items))
    
    #print("brain first 10 keys:", list(grn)[:10])
    degset = {row[0] for row in degs} # location of DEG col in input chart
    #print("deg first 10 keys:", list(degset)[:10])
    print("overlap:", len(set(grn) & degset))

    # merge list of lists with brain grn, cutting the grn down to just the keys in disorders
    """
    1. deg_grn_tfsonly preserves network topology and adds annotation (has detf regulating all genes). 
    A transcription factor's expression is abnormal in a disorder, so its regulatory output may be altered. 
    A TF that is upregulated or downregulated can shift expression of many genes that themselves are not significantly DE.
    Common where a master regulator shifts subtly but downstream genes do not pass DEG thresholds.
    Detects: regulatory amplification, TFs with unusually many DEG targets, TFs whose downstream genes shift collectively but weakly
    fraction DEG per TF = (# DEG targets) / (# total targets): measures how strongly that TF's regulated genes overlaps with disease signal.
    
    2. deg_grn_both performs two-stage filtering, producing a smaller induced subgraph over DEG nodes (detf regulates deg)
    Thus, identifying direct dysregulated regulatory cascades eg synaptic regulation, immune activation, chromatin modifiers
    Gives small but highest-confidence regulatory modules useful for pathway discovery, disease modules, candidate causal regulators
    Limitation is that many regulatory effects are lost because DEG thresholds are harsh, 
    brain tissue averages many cell types,
    TF activity does not always correlate with TF expression
    
    3. deg_grn_genesonly has tf regulating degs. TF activity changes without expression change.
    post-translational TF activation (phosphorylation), chromatin accessibility change, 
    cofactor-dependent TF activation, 
    cell-type composition changes    
    """
    
    detfdeggrn = gu.de_grn_both(grn,degs) # produces grn of differential tfs AND differential gene targets
    tfgrn = gu.de_grn_tfsonly(grn,degs) # produces grn of differential tfs only
    deggrn = gu.de_grn_degsonly(grn,degs) # produces a grn of differential gene targets only
    
    # turn into edgelist
    detfdeggrnedgelist = gu.adjlist2edgelist(detfdeggrn)
    print("size of deg-detf-grn: ", len(detfdeggrn)) # graph this somehow comparing all of the above?
    print("size of tf-grn: ", len(tfgrn)) 
    print("size of deg-grn: ", len(deggrn)) 
    
    """for tf, edges in deggrn.items(): 
        for e in edges[:5]:
            print(len(e))
        break
    print("deg-grn keys:", list(deggrn)[:10])
    """
    
    # run stats from graph_algos here!
    
    # check number of positive vs negative DETFs
    ga.edgeweight_summary(brains)
    ga.edgeweight_summary(grn)
    ga.edgeweight_summary(brainother)
    ga.edgeweight_summary(brainbg)
    ga.edgeweight_summary(detfdeggrn)
    ga.edgeweight_summary(tfgrn)
    
    ga.log2fc_summary(detfdeggrn)
    ga.log2fc_summary(tfgrn)
    
    """
    A TF ranks highly when it regulates many DEGs, the edges have strong PANDA weights and those genes show large expression changes
    = candidates for regulatory drivers of the disorder transcriptome.
    
    TFs with high normalized score regulate dense clusters of dysregulated genes.
    """
    bddrivers = ga.regulator_detection(grn, gu.disorder_list('data/DEGDataSample.csv','BD'))
    szdrivers = ga.regulator_detection(grn, gu.disorder_list('data/DEGDataSample.csv','SZ'))
    bdszdrivers = ga.regulator_detection(grn, gu.disorder_list('data/DEGDataSample.csv','SZ','BD'))
    
    for r in bddrivers[:10]:
        print(r,"\n")
    
    

    """G = nx.from_dict_of_lists(
    {k: [x[0] for x in v] for k, v in deggrn.items()},
    create_using=nx.DiGraph())"""
    #nx.readwrite.json_graph.cytoscape_data(G)
    #pyc.create_network_from_networkx(G, title="Test DEGs")
    
    
    
    # visualize graph in new file
    gv.viz_graph(deggrnedgelist,'results/deggrn.html')
    gv.visualize_deg_grn(detfdeggrn)
    gv.pyviz_deggrn(detfdeggrn)
    
    ## later on, add variable that tracks the names of the disorders being indexed so I can create a file with that name
                      
    # calculate graph metrics (use an existing package?)
    
    return 

if __name__ == '__main__':
    main()

