import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from bct import motif3struct_bin, motif3funct_bin, motif4struct_bin, motif4funct_bin
from connectivity_graph import compute_adjacency, load_matrix, load_conn_graph, p1_5

PLOTS = False

""" TO INSTALL bctpy 
$ git clone https://github.com/aestrivex/bctpy
$ cd bctpy
$ python3 setup.py build
$ python3 setup.py install
$ sudo mv motif34lib.mat /usr/local/lib/python3.6/dist-packages/bctpy-0.5.2-py3.6.egg/bct
"""


def triadSignificanceProfile(G, triad_cfg="fuffa"):
    """
    Compute the significance profile of the patterns mapped in triad_cfg, 
    inside directed graph G.
        - G          : directed graph representing the network; 
        - triads_cfg : dict mapping interesting triadic patterns codes, 
            as in nx.triadic_census(), with explicit names. 
            (e.g. triad_cfg = {'003' : 'Null', '012' : 'Single-edge'})
    """
    M = compute_adjacency(load_matrix())
    m_3, M_3 = motif3funct_bin(M)

    in_degree_sequence  = [d for n, d in G.in_degree()]   # in degree sequence
    out_degree_sequence = [d for n, d in G.out_degree()]  # out degree sequence
    #print("In_degree sequence %s" % in_degree_sequence)
    #print("Out_degree sequence %s" % out_degree_sequence)

    random_nets_census = []
    for i in range(100):
        rand_G = nx.directed_configuration_model(in_degree_sequence, out_degree_sequence, create_using=nx.DiGraph, seed=i)
        adj_M  = nx.to_numpy_array(rand_G)

        rand_m_3, rand_M_3 = motif3funct_bin(adj_M)
        random_nets_census.append(rand_m_3)

    real_census   = m_3
    random_census = random_nets_census
    print(real_census)
    print(random_census)

    z_score = []
    for p in range(len(real_census)):
        print(p)
        N_real_p = real_census[p]
        N_rand_p = np.mean(random_census[p])
        std = np.std(random_census[p])

        z_p =  ((N_real_p - N_rand_p)/std if std != 0 else 0)
        z_score.append(z_p)

    sp = []
    for i in range(len(z_score)):
        z_norm = np.linalg.norm(z_score)
        norm_z_score = (z_score[i]/z_norm if z_norm != 0 else z_score[i])
        sp.append(round(norm_z_score,4))

    print(sp)
    return sp


def p3_1():
    M = compute_adjacency(load_matrix())
    m_3, M_3 = motif3funct_bin(M)

    print("[3.1] >> Motif frequency:", m_3)
    plt.bar(np.arange(1, 14), m_3)
    plt.xlabel("Motif ID")
    plt.ylabel("frequency")
    plt.xticks(np.arange(0, 14, 1))
    plt.title("Network class-3 motif frequency in the graph")
    plt.show()

    print("[3.1] >> Motif 1 node frequency:", M_3[0])
    plt.matshow(M_3)
    plt.xlabel("Node ID")
    plt.ylabel("Motif ID")
    plt.title("Node class-3 motif frequency fingerprint")
    plt.show()
    return M_3


def p3_2(G):
    """
    Plots a new graph with the same nodes as G, containing only G edges involved
    in motif of type 1.
    """
    print("[3.2] >> found {} edges between {} nodes in the original graph.".format(len(G.edges()), len(G.nodes())))
    motif_G = nx.create_empty_copy(G)
    for node in G.nodes():
        for e1 in G.in_edges(node):
            for e2 in G.in_edges(node):
                if e2 != e1 and (e1[0], e2[0]) not in G.edges() and (e2[0], e1[0]) not in G.edges() and (e1[1], e1[0]) not in G.edges() and (e2[1], e2[0]) not in G.edges():
                    motif_G.add_edge(e1[0],e1[1], color='r')
                    motif_G.add_edge(e2[0],e2[1], color='r')
    for edge in G.edges():
        if edge not in motif_G.edges():
            motif_G.add_edge(edge[0], edge[1], color='b')
    print("[3.2] >> found {} edges between {} nodes in the new graph.".format(len(motif_G.edges()),len(motif_G.nodes())))
    p1_5(motif_G, point='3.2')


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


def p3_4():
    """
    Plots frequency of class-4 motif involving in the graph 
    described by 'adj_mat' adjacency matrix.
    """
    M = compute_adjacency(load_matrix())
    m_4, M_4 = motif4funct_bin(M)
    print("[3.4] >> Displaying frequencies of 4-node motifs.")
    plt.bar(np.arange(0,199), m_4)
    plt.subplots_adjust(left=0.05, right=0.95, top=0.92, bottom=0.08)
    plt.xlabel("Motif ID")
    plt.ylabel("frequency")
    ticks = np.arange(0, 200, 5)
    ticks[0] = 1
    ticks[-1] = 199
    plt.xticks(ticks, rotation='vertical')

    """
    makes the plot larger to fit all 200 values, but it's unpractical

        plt.gca().margins(x=0)
        plt.gcf().canvas.draw()
        maxsize = 11
        m=0.2
        s = maxsize/plt.gcf().dpi*200+2*m
        margin = m/plt.gcf().get_size_inches()[0]
        plt.gcf().subplots_adjust(left=margin, right=1.-margin)
        plt.gcf().set_size_inches(s, plt.gcf().get_size_inches()[1])
    """
    
    plt.title("Network class-4 motif frequency in the graph")
    plt.show()

    return m_4, M_4


G = load_conn_graph(conn="pdc", freq=10, run="R01")
#p3_2(G)

sp = triadSignificanceProfile(G)
print(G)
