import matplotlib.pyplot as plt
import numpy as np
import algos
# by iris

## code below is referenced from either previous assignments or matplotlib visualization references
# for both box and violin, clustcoeff, closeness, eccentricity not showing up fully

def net_compare(m_adjlist, e_adjlist, outfile1, outfile2):
    
    # call all functions from thomas's algos code
    ## mesenchymal
    m_cluscoeff = algos.calculate_clustering_coeff(algos.unweight_adjlist(m_adjlist, algos.get_nodes(m_adjlist))) 
    m_closeness = algos.closeness_centralities(m_adjlist)
    m_eccentricity = algos.eccentricity_centralities(m_adjlist)
    m_betweenness = algos.betweenness_centralities(m_adjlist, True)
    
    ## epithelial
    e_cluscoeff = algos.calculate_clustering_coeff(algos.unweight_adjlist(e_adjlist, algos.get_nodes(e_adjlist))) 
    e_closeness = algos.closeness_centralities(e_adjlist)
    e_eccentricity = algos.eccentricity_centralities(e_adjlist)
    e_betweenness = algos.betweenness_centralities(e_adjlist, True)
    
    
    # list of dictionaries of lists initializing the network measures for epithelial vs mesenchymal
    measures = {
        'Epithelial': {'degree': [], 'clustering': [], 'closeness': [], 'eccentricity': [], 'betweenness': []}, 
        'Mesenchymal': {'degree': [], 'clustering': [], 'closeness': [], 'eccentricity': [], 'betweenness': []}
    }

    for node in m_adjlist:
        measures['Mesenchymal']['degree'].append(algos.weighted_node_degree(algos.reciprocals(m_adjlist), node)) ##setting the values of each net. measure per node status
        measures['Mesenchymal']['clustering'].append(m_cluscoeff[node])
        measures['Mesenchymal']['closeness'].append(m_closeness[node])
        measures['Mesenchymal']['eccentricity'].append(m_eccentricity[node])
        measures['Mesenchymal']['betweenness'].append(m_betweenness[node])

    for node in e_adjlist:
        measures['Epithelial']['degree'].append(algos.weighted_node_degree(algos.reciprocals(e_adjlist), node)) ##setting the values of each net. measure per node status
        measures['Epithelial']['clustering'].append(e_cluscoeff[node])
        measures['Epithelial']['closeness'].append(e_closeness[node])
        measures['Epithelial']['eccentricity'].append(e_eccentricity[node])
        measures['Epithelial']['betweenness'].append(e_betweenness[node])
    ## faceted boxplots to compare node measures 
    """BOX PLOTS"""

    ##separating data into graph measures
    deg_data = [measures['Epithelial']['degree'], measures['Mesenchymal']['degree']] 
    clust_data = [measures['Epithelial']['clustering'], measures['Mesenchymal']['clustering']] ##need to check if this is how the categories work for the cell data
    close_data = [measures['Epithelial']['closeness'], measures['Mesenchymal']['closeness']]
    eccen_data = [measures['Epithelial']['eccentricity'], measures['Mesenchymal']['eccentricity']]
    betw_data = [measures['Epithelial']['betweenness'], measures['Mesenchymal']['betweenness']]

    fig1, axs = plt.subplots(1, 5, layout = 'constrained')

    axs[0].boxplot(deg_data, labels=['E', 'M'],
        patch_artist=True, boxprops={'facecolor': 'lightcoral'})
    axs[0].set_title('Degree') 
    
    axs[1].boxplot(clust_data, labels=['E', 'M'],
        patch_artist=True, boxprops={'facecolor': 'gold'})
    axs[1].set_title('Clust. Coefficient')

    axs[2].boxplot(close_data, labels=['E', 'M'],
        patch_artist=True, boxprops={'facecolor': 'yellowgreen'})
    axs[2].set_title('Closeness')

    axs[3].boxplot(eccen_data, labels=['E', 'M'],
        patch_artist=True, boxprops={'facecolor': 'skyblue'})
    axs[3].set_title('Eccentricity')

    axs[4].boxplot(betw_data, labels=['E', 'M'],
        patch_artist=True, boxprops={'facecolor': 'mediumpurple'})
    axs[4].set_title('Betweenness')

    fig1.suptitle('Boxplots of Connectivity Measures for E vs. M Cells') 

    fig1.savefig(outfile1) 
    plt.close(fig1)

    """VIOLIN PLOTS"""
    ## faceted violin plots for showing distribution of individual points compared to boxplots below
    fig2, axs = plt.subplots(1, 5, layout = 'constrained') # are there other options other than constrained?

    # plot violin plots
    axs[0].violinplot(deg_data, showmeans=False, showmedians=True)
    axs[0].set_title('Degree')

    axs[1].violinplot(clust_data, showmeans=False, showmedians=True)
    axs[1].set_title('Clust. Coefficient')

    axs[2].violinplot(close_data, showmeans=False, showmedians=True)
    axs[2].set_title('Closeness')

    axs[3].violinplot(eccen_data, showmeans=False, showmedians=True)
    axs[3].set_title('Eccentricity')

    axs[4].violinplot(betw_data, showmeans=False, showmedians=True)
    axs[4].set_title('Betweenness')


    
    for ax in axs:
        tk = [1,2]
        ax.set_xticks(ticks=tk,labels=['E', 'M'])
    
    
    # ax = axs[0]
    # ax.yaxis.grid(True)
    # tk = []
    # while len(deg_data) != len(tk): tk.append(1)
    # print(tk)
    # ax.set_xticks(ticks=tk,labels=['Epithelial', 'Mesenchymal'])
    

    fig2.suptitle('Violin Plots of Connectivity Measures for E vs. M Cells')

    fig2.savefig(outfile2) 
    plt.close(fig2)


    #return m_cluscoeff, m_closeness, m_eccentricity, m_betweenness, deg_data
    return deg_data, clust_data, close_data, eccen_data, betw_data








