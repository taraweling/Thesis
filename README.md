# Cross-Disorder GRN
For my thesis, I use the GRAND database which implements PANDA (https://netzoo.github.io/zooanimals/panda/) to produce tissue-specific bipartite Gene-Regulatory Networks (GRNs).

I subsequently cross-reference DEGs for various disorders to TFs in the GRN to produce relevant pathways. 

Code is included in src. The main function is in run.py, which calls other functions in that file, while helper functions that are used repeatedly are in utils.py.

Since the processes rely on resource-intensive, time-consuming loops, they are optimally run on computing nodes with the settings found in ther sbatch files.
