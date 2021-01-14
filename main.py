from connectivity_graph import *
from motif_analysis import *
from community_detection import *

G = load_conn_graph(conn="pdc", freq=10, run="R01")


# PART 1
p1_1()

p1_2()
p1_3()
p1_4()

if PLOTS:
    print("[1.5] >> {} nodes\t{} edges".format(len(G.nodes()), len(G.edges())))
    p1_5(G)



# PART 2




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
