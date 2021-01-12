import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from bct import motif3struct_bin, motif4struct_bin
from connectivity_graph import compute_adjacency, load_matrix, load_conn_graph, p1_5

PLOTS = False

""" TO INSTALL bctpy 
$ git clone https://github.com/aestrivex/bctpy
$ cd bctpy
$ python3 setup.py build
$ python3 setup.py install
$ sudo mv motif34lib.mat /usr/local/lib/python3.6/dist-packages/bctpy-0.5.2-py3.6.egg/bct
"""

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
    print("[3.2] >> found {} edges between {} nodes.".format(len(G.edges()), len(G.nodes())))
    p1_5(motif_G)

def p3_3(freq_mat, ch_name="Po4"):
    """
    Plots frequency of motif involving 'ch_name' channel.
    """
    # TODO ch. 59 <-> Po4.
    print("[3.3] >> displaying the motif in which {} is involved.".format(ch_name))
    plt.bar(np.arange(1, 14), freq_mat[:, 59])
    plt.xlabel("Motif ID")
    plt.ylabel("frequency")
    plt.xticks(np.arange(0,14,1))
    plt.title("Motif frequency - channel {}".format(ch_name))
    plt.show()

def p3_4(adj_mat):
    """
    Plots frequency of class-4 motif involving in the graph 
    described by 'adj_mat' adjacency matrix.
    """
    m_4, M_4 = motif4struct_bin(adj_mat)

    plt.bar(np.arange(1,200), m_4)
    plt.xlabel("Motif ID")
    plt.ylabel("frequency")
    plt.xticks(np.arange(0,200,1))
    plt.title("Network class-4 motif frequency in the graph")
    plt.show()

    return m_4, M_4



M = adj_mat = compute_adjacency(load_matrix())
m_3, M_3 = motif3struct_bin(M)

if PLOTS:
    print("Motif frequency:",m_3)
    plt.bar(np.arange(1,14), m_3)
    plt.xlabel("Motif ID")
    plt.ylabel("frequency")
    plt.xticks(np.arange(0,14,1))
    plt.title("Network class-3 motif frequency in the graph")
    plt.show()

    print("Motif 1 node frequency:",M_3[0])
    plt.matshow(M_3)
    plt.xlabel("Node ID")
    plt.ylabel("Motif ID")
    plt.title("Node class-3 motif frequency fingerprint")
    plt.show()



#### 3.2
p3_2(G=load_conn_graph())

#### 3.3
p3_3(freq_mat=M_3)

#### 3.4
#p3_4(adj_mat=adj_mat)
