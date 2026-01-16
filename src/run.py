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
    testadjlist = utils.make_adjlist('data/test.csv',0) # truong paper uses threshold of 0.7
    print(testadjlist)
    print('\n')
    
    # get the GRAND GRN by running on miles' nodes
    #brainother = utils.make_adjlist('data/Brain_Other.csv',0.0) 
    #brainbg = utils.make_adjlist('data/Brain_Basal_Ganglia.csv',0.0)
    
    # test merge_adjlist's weight-averaging feature: AHR-A1BG SHOULD have a weight of ~0.2
    testadjlist2 = {'AHR': [('A1BG', 0.0), ('A4GALT', 0.314864), ('AHR', 0.131604), ('ALX3', 0.208409)]}
    testmerge = utils.merge_adjlist(testadjlist,testadjlist2)
    #print(testmerge,'\n')

    # merge the two brain GRNs (also on cluster)
    #brain = utils.merge_adjlist(brainother, brainbg)
    
    # get adjlists per disorder 
    #degs_disorder outputs list
    testdisorder = ['AHR', 'ALX3']
    
    # reduce GRN to just those in disorder list
    deggrn = utils.filter_adjlist(testmerge,testdisorder)
  
    # convert to ENSEMBL
    ensemblify(deggrn)
                

    
    print(deggrn)
    
    ## for now, input is DEGDataSample.csv
    
    # visualize graph in new file
    
    
    # calculate graph metrics (use an existing package?)
    
    
    return 



def disorder_lists(*disorder:str): # generates the concatenated list of n number of disorders queried by terminal input
    
    
    
    
    return 


def ensemblify(deggrn):
    for k in list(deggrn): # loops through each tf and its values
        new = utils.gene2ensembl(k)
        deggrn[new] = deggrn.pop(k) # ensemblfies all keys
        v = deggrn[new] # looks at each TF's targets
        for i, gene in enumerate(v):
            newgene = utils.gene2ensembl(gene[0])
            v[i] = (newgene,gene[1])
    return deggrn

if __name__ == '__main__':
    main()

