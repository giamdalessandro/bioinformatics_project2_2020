from commons import *
from graph_indices import *
from motif_analysis import *
from connectivity_graph import *
from community_detection import *

"""
TODO
- 1.4   try to find the libraries
"""

# we check if we can use just one frequencies of the alpha band or the entire band
check_mean_var_EEG_rithm(plot=False)

# since the values of the variance matrix are two order of magnitude less than the mean matrix, 
# we will use a specific frequence - 10Hz
freq = 10
conn = 'pdc'

####### TASK 1 ########


# NOTE that we compute the matrix one time and then we save them
# to re-run the code that computes them, delete the matrix files
# e.g. data/pdc_R01_10Hz_auto.txt
# e.g. data/dtf_R02_12Hz_auto.txt

for r in RUNS:
    print("\n============================== TASTK 1 - {} ==============================".format(r))
    # p1_1(file_name="data/S003{}_fixed".format(r), freq=freq, run=r)  #NOTE: p1_1 performs also point 1.2
    # p1_3(conn_method=conn, freq=10, run=r)
    # p1_4()
    # p1_5(load_conn_graph(conn=conn, freq=freq, run=r))
    # p1_6(file_name="data/S003{}_fixed".format(r), freq=25, run=r)

# a

####### TASK 2 ########
#cf, pl = cf_pl()       # speeds up p2_4
for r in RUNS:
    print("\n============================== TASTK 2 {} ==============================".format(r))
    #cf_real, pl_real = p2_1(conn=conn, freq=freq, run=r)
    #small_worldness  = p2_2(cf_real, pl_real)
    #p2_3(freq=freq, run=r)
    #p2_4(run=r, cf_rand=cf, pl_rand=pl)
    #p2_5(load_conn_graph(conn=conn, freq=freq, run=r))
    p2_6(run=r)
    #p2_7(run=r)



####### TASK 3 ########
'''
for r in RUNS:
    print("\n============================== TASK 3 - {} ==============================".format(r))
    M_3 = p3_1()
    p3_2(G)
    p3_3(freq_mat=M_3)
    p3_4()
'''


####### TASK 4 ########
'''
for r in RUNS:
    print("\n============================== TASK 4 - {} ==============================".format(r))
    partition, partition_louvain = p4_1(run=r)
    p4_2(run=r, partition=partition)
    partition_infomap = p4_3(run=r)

    for S1 in partition_louvain.values():
        for S2 in partition_infomap:
            j = jaccard(S1, S2)             # S1 and S2 are list of nodes
            if j > 0:
                print("Louvain community:", S1)
                print("Infomap community:", S2)
                print("Jaccard is {}%\n".format(j))

'''


# to see the mapping between integers and channels in latex
# mapping = map_index_to_channels()
# for k,v in mapping.items():
#     print("{} & {} \\\\ \hline".format(k,v))