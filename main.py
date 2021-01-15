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
conn_method = 'pdc'

####### TASK 1 ########


# NOTE that we compute the matrix one time and then we save them
# to re-run the code that computes them, delete the matrix files
# e.g. data/pdc_R01_10Hz_auto.txt
# e.g. data/dtf_R02_12Hz_auto.txt

p1_1(file_name="data/S003R01_fixed", freq=freq, run='R01')  # NOTE: p1_1 performs also point 1.2
p1_3(conn_method=conn_method, freq=10, run='R01')
p1_4()
p1_5(load_conn_graph(conn=conn_method, freq=freq, run="R01"))
p1_6(file_name="data/S003R01_fixed", freq=25, run='R01')

p1_1(file_name="data/S003R02_fixed", freq=freq, run='R02')  # NOTE: p1_1 performs also point 1.2
p1_3(conn_method=conn_method, freq=10, run='R02')
p1_4()
p1_5(load_conn_graph(conn=conn_method, freq=freq, run="R02"))
p1_6(file_name="data/S003R02_fixed", freq=25, run='R02')


####### TASK 2 ########

#N = 64
#nodes_idx = [i for i in range(N)]  # [0,1,2..., N]
#max_edges = N*(N-1)

# get adjacency matrix
conn_mat = load_matrix(conn_method=conn_method)
print("mat shape:", conn_mat.shape)

threshold_20_percent_density = 0.137 # 0.04597 for PDC
adj_mat = compute_adjacency(conn_mat, threshold=threshold_20_percent_density)  
print("Resulting network density:", np.sum(adj_mat)/4032)

print('\n================== P 2.1 ==============================')
Cf_real, PL_real = graph_indices_part_2_1(adj_mat)


print('\n================== P 2.2 ==============================')
print('small worls formula = (Cf_G/Cf_rand)/(PL_G/PL_rand)')
# small worls formula = (Cf_G/Cf_rand)/(PL_G/PL_rand)
Small_worldness = graph_indices_part_2_2(Cf_real, PL_real, random_graph='erdos')    # or 'watts'

print('\n================== P 2.3 ==============================')
p2_3(freq=10, run='R01', threshold_pdc=0.1226, threshold_dtf=0.1378)


print('\n================== P 2.4 ==============================')
# just change threshold values in crearting adjacency matrix to tune density as 
# mentioned in P 1.3
#  densities = [1%, 5%, 10%, 20%, 30%, 50%]
thresholds = [0.41, 0.24, 0.187, 0.137, 0.1, 0.055]
graph_indices_part_2_4(conn_mat, thresholds)

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
