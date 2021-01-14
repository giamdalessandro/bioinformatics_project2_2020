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

def task1():
    # we check if we can use just one frequencies of the alpha band or the entire band
    check_mean_var_EEG_rithm()

    # since the values of the variance matrix above thresold are few, we will use a specific frequence - 10Hz
    p1_1()
    p1_2()
    p1_3()
    p1_4()
    p1_5(G)

 
G = load_conn_graph(conn="pdc", freq=10, run="R01")
print("[1.5] >> {} nodes\t{} edges".format(len(G.nodes()), len(G.edges())))

task1()

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
