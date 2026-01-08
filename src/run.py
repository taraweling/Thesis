import csv
from pyvis.network import Network as net
import networkx as nx
import matplotlib.pyplot as plt
import utils

def main():
    
    # obtain adjlist for each disorder (handle in main or helper fn?)
    ## what information do I want in my graph?
    """
    as a bipartite graph, node attribute should either be TF or gene. 
    the TFs are DEGs and the target genes can be either DEGs or not-differentially expressed genes in the GRAND regulatory network 
    """
    
    # test make_adjlist
    testadjlist = utils.make_adjlist('test.csv',0) # truong paper uses threshold of 0.7
    print(testadjlist)
    print('\n')
    # get the GRAND GRN by running on miles' nodes
    #brainother = utils.make_adjlist('Brain_Other.csv',0.0) #
    #brainbg = utils.make_adjlist('Brain_Basal_Ganglia.csv',0.0)
    
    # test merge_adjlist's weight-averaging feature: AHR-A1BG SHOULD have a weight of ~0.2
    newadjlist = {'AHR': [('A1BG', 0.0), ('A4GALT', 0.314864), ('AHR', 0.131604), ('ALX3', 0.208409)], 'AIRE': [('A1CF', 3.340785), ('A2M', 0.365588)], 'ALX1': [('A1CF', 0.181709), ('A2M', 0.485148)], 'ALX3': [('A1CF', 0.208687), ('A2M', 3.764654)]}
    testmerge = utils.merge_adjlist(testadjlist,newadjlist)
    print(testmerge)

    # merge the two brain GRNs (also on cluster)
    #brain = utils.merge_adjlist(brainother, brainbg)
    
    # get adjlists per disorder 
    testlist = ['AHR', 'ALX3','A1BG']
    ## reduce brain 
    # A1BG shouldn't appear in
    
    ## for now, input is DEGDataSample.csv
    
    
    # function to combine adjlists?
    
    # create graph
    
    # calculate graph metrics (use an existing package)
    
    #                                                                                                                  
    # visualize
    
    
    return 

def disorder_lists(*disorder:str): # generates the concatenated list of n number of disorders 
    
    
    
    
    return 

if __name__ == '__main__':
    main()

