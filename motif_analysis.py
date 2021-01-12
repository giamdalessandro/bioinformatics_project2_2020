import numpy as np
import networkx as nx
from connectivity_graph import load_conn_graph

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
plt.bar(np.arange(1,14), m)
plt.xlabel("Motif ID")
plt.ylabel("frequency")
plt.xticks(np.arange(0,14,1))
plt.title("Network motif frequency in the graph")
plt.show()

# che cazz è
print(m)
plt.matshow(M)
plt.xlabel("Node ID")
plt.ylabel("Motif ID")
plt.title("Node motif frequency fingerprint")
plt.show()

#### 3.2
G = load_conn_graph()
empty_G = nx.create_empty_copy(load_conn_graph())

