from commons import *
from graph_indices import *
from motif_analysis import *
from connectivity_graph import *
from community_detection import *

"""
TODO
- 1.1   shuffle rows/cols?
- 1.4   try to find the libraries
- 3.1   check that we are performing the statistical test
"""

# we check if we can use just one frequencies of the alpha band or the entire band
check_mean_var_EEG_rithm(plot=False)

# since the values of the variance matrix above thresold are few, we will use a specific frequence - 10Hz
freq = 10
conn = 'pdc'

####### TASK 1 ########


# NOTE that we compute the matrix one time and then we save them
# to re-run the code that computes them, delete the matrix files
# e.g. data/pdc_R01_10Hz_auto.txt
# e.g. data/dtf_R02_12Hz_auto.txt
'''
for r in RUNS:
    print("\n============================== {} ==============================".format(r))
    print("'{}'".format(r))
    p1_1(file_name="data/S003{}_fixed".format(r), freq=freq, run=r)  # NOTE: p1_1 performs also point 1.2
    p1_3(conn_method=conn, freq=10, run=r)
    p1_4()
    p1_5(load_conn_graph(conn=conn, freq=freq, run=r))
    p1_6(file_name="data/S003{}_fixed".format(r), freq=25, run=r)
'''


####### TASK 2 ########

for r in RUNS:
    print("\n============================== {} ==============================".format(r))
    cf_real, pl_real = p2_1(conn=conn, freq=freq, run=r)
    small_worldness  = p2_2(cf_real, pl_real, random_graph='erdos')  # or 'watts'
    p2_3(freq=freq, run=r)
    p2_4(run=r)

'''


print('\n================== P 2.5 ==============================')
G = load_conn_graph(conn="pdc", freq=10, run="R01")
p2_5(G)

print('\n================== P 2.6 ==============================')
p2_6('R01')
p2_6('R02')

print('\n================== P 2.7 ==============================')
threshold = threshold_20_percent_density
graph_indices_part_2_7(conn_mat, threshold)

'''




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
