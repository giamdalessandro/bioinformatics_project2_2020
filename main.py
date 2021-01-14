from connectivity_graph import *
from motif_analysis import *
from community_detection import *

"""
TODO
- 1.1   shuffle rows/cols?
- 1.3
- 1.4   try to find the libraries
- 1.6   same as 1.1
- 2.6
- 3.1   struct or func motif? 
        check that we are performing the statistical test
"""

# we check if we can use just one frequencies of the alpha band or the entire band
check_mean_var_EEG_rithm(plot=False)

# since the values of the variance matrix above thresold are few, we will use a specific frequence - 10Hz
freq = 10


####### TASK 1 ########


# NOTE that we compute the matrix one time and then we save them
# to re-run the code that computes them, delete the matrix files
# e.g. data/pdc_R01_10Hz_auto.txt
# e.g. data/dtf_R02_12Hz_auto.txt

p1_1(file_name="data/S003R01_fixed", freq=freq, run='R01')  # NOTE: p1_1 performs also point 1.2
p1_3(conn_method='pdc', freq=10, run='R01')




p1_1(file_name="data/S003R02_fixed", freq=freq, run='R02')  # NOTE: p1_1 performs also point 1.2
p1_3(conn_method='pdc', freq=10, run='R02')




#G = load_conn_graph(conn="pdc", freq=freq, run="R01")
#print("Loaded {} nodes\t{} edges".format(len(G.nodes()), len(G.edges())))


'''
# PART 3
M_3 = p3_1()
p3_2(G)
p3_3(freq_mat=M_3)
p3_4()



# PART 4
partition, partition_louvain = p4_1()
p4_2(G, partition)
partition_infomap = p4_3(G)

#NOTE: implementare una metrica (Jaccard ?) per vedere quanto le due partition sono simili

for S1 in partition_louvain.values():
    for S2 in partition_infomap:
        j = jaccard(S1, S2)             # S1 and S2 are list of nodes
        if j > 0:
            print("Louvain community:", S1)
            print("Infomap community:", S2)
            print("Jaccard is {}%\n".format(j))
'''
