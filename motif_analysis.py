import numpy as np
import networkx as nx
from connectivity_graph import load_conn_graph, p1_5

"""
triad_cfg = {
	'021D': 'type-1',
	'021C': 'Three-chain',
	'021U': 'type-3',
    '111D': 'type-4',
    '111U': 'type-5',
    '030T': 'Feed-forward',
	'030C': 'Feedback',
    '201' : 'type-8',
    '120D': 'type-9',
    '120U': 'type-10',
    '120C': 'type-11',
    '210' : 'type-12',
    '300' : 'type-13'
}

net_G = load_conn_graph()
census = nx.triadic_census(net_G)

f_census = {}
print('\ncomputing network triadic census...\n')
print('triadType  \tN')
print('--------------------')
for k,v in sorted(census.items()):
	if k in triad_cfg:
		f_census[triad_cfg[k]] = [v]
		print(triad_cfg[k] + ': \t' + str(v))
"""



"""
$ git clone https://github.com/aestrivex/bctpy
$ cd bctpy
$ python3 setup.py build
$ python3 setup.py install
$ sudo mv motif34lib.mat /usr/local/lib/python3.6/dist-packages/bctpy-0.5.2-py3.6.egg/bct
"""

import matplotlib.pyplot as plt
from bct import motif3struct_bin
from connectivity_graph import compute_adjacency, load_matrix

M = adj_mat = compute_adjacency(load_matrix())
m, M = motif3struct_bin(M)
print("Motif frequency:",m)
plt.bar(np.arange(1,14), m)
plt.xlabel("Motif ID")
plt.ylabel("frequency")
plt.xticks(np.arange(0,14,1))
plt.title("Network motif frequency in the graph")
plt.show()

# che cazz è
print("Motif 1 node frequency:",M[0])
plt.matshow(M)
plt.xlabel("Node ID")
plt.ylabel("Motif ID")
plt.title("Node motif frequency fingerprint")
plt.show()



#### 3.2
def p3_2(G):
    """
    Plots a new graph with the same nodes as G, containing only G edges involved
    in motif of type 1.
    """
    motif_G = nx.create_empty_copy(G)

    for node in G.nodes():
        for e1 in G.in_edges(node):
            for e2 in G.in_edges(node):
                if e2 != e1 and (e1[0],e2[0]) not in G.edges() and (e2[0],e1[0]) not in G.edges():
                    motif_G.add_edge(e1[0],e1[1])
                    motif_G.add_edge(e2[0],e1[1])

    p1_5(motif_G)

#### 3.3
def p3_3(freq_mat, ch_name="Po4"):
    # ch. 59 <-> Po4.
    plt.bar(np.arange(1,14), M[:,59])
    plt.xlabel("Motif ID")
    plt.ylabel("frequency")
    plt.xticks(np.arange(0,14,1))
    plt.title("Motif frequency - channel {}".format(ch_name))
    plt.show()



#p3_2(G=load_conn_graph())
p3_3(freq_mat=M)