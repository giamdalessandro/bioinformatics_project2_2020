from connectivity_graph import *
from motif_analysis import *

# PART 1

p1_1()

G = load_conn_graph(conn="pdc", freq=10, run="R01")

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
