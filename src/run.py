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

    # merge the two brain GRNs (also on cluster)
    #brainother = gu.ensemblify(brainother)
    #brainbg = gu.ensemblify(brainbg)
    brains= gu.merge_adjlist(brainother, brainbg)
    grn = gu.ensemblify(brains)
    
    """Get the first key
    first_node = next(iter(brains))
    print("non ensemble",first_node)"""

    # should I do a bunch of stats calls on the merged brain grn to compare?
    
    
    # get adjlists per disorder by inputting the name of the disorder CURRENTLY ENSEMBLIFY IS TAKING TOO LONG
    ## (options = AD, ADHD, BD, SZ, MDD, OCD) and the file location, outputting a list of lists   
    test = gu.disorder_list('data/DEGDataSample.csv','BD','SZ')
    
    # lambda function that creates a list from the keys in the disorders dict?
    
    # convert geneIDs in disorders list of lists to ENSEMBL IDs
    degs = gu.ensemblifylist(test)
    
    # reduce GRN to just keys from disorders 
    #deggrn = gu.filter_adjlist(brain,disorders)
    
    # merge list of lists with brain grn, cutting the grn down to just the keys in disorders
    print(len(grn),len(degs))
    degset = {row[4] for row in degs} # location of DEG col

    print("brain keys:", list(grn)[:1000])
    print("degset:", list(degset)[:10])
    print("overlap:", len(set(grn) & degset))

    deggrn = gu.deg_grn_tfsonly(grn,degs)
    print("size of deg-grn: ", len(deggrn))
    # run stats from graph_algos here!
    
    """G = nx.from_dict_of_lists(
    {k: [x[0] for x in v] for k, v in deggrn.items()},
    create_using=nx.DiGraph())"""
    #nx.readwrite.json_graph.cytoscape_data(G)
    #pyc.create_network_from_networkx(G, title="Test DEGs")
    
    # turn into edgelist
    deggrnedgelist = gu.adjlist2edgelist(deggrn)
    
    # visualize graph in new file
    gv.viz_graph(deggrnedgelist,'results/disorder.html')
    gv.visualize_deg_grn(gu.deg_grn_both(grn,degs))
    #gv.lastditch(deggrn)
    ## later on, add variable that tracks the names of the disorders being indexed so I can create a file with that name
                      
    # calculate graph metrics (use an existing package?)
    
    return 


if __name__ == '__main__':
    main()

