import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from bct import motif3struct_bin, motif4struct_bin
from connectivity_graph import compute_adjacency, load_matrix, load_conn_graph, p1_5

PLOTS = False


def significanceProfile(M, nrand=100):
    """
    Compute the significance profile of the motifs in a directed graph G, represented by the adjacency matrix M.
        - M : adjacency matrix representing the network; 
    """
    m_3, M_3 = motif3struct_bin(M)

    in_degree_sequence  = [d for n, d in G.in_degree()]   # in-degree sequence
    out_degree_sequence = [d for n, d in G.out_degree()]  # out-degree sequence

    random_nets_census = []
    over_rep = []
    for i in range(nrand):
        rand_G = nx.directed_configuration_model(in_degree_sequence, out_degree_sequence, create_using=nx.DiGraph, seed=i)
        adj_M  = nx.to_numpy_array(rand_G)

        rand_m_3, rand_M_3 = motif3struct_bin(adj_M)
        random_nets_census.append(rand_m_3)
        ov = []
        for i in range(len(rand_m_3)):
            if rand_m_3[i] > m_3[i]: 
                ov.append(1)
            else:
                ov.append(0)

        over_rep.append(ov)

    real_census   = m_3
    random_census = np.array(random_nets_census)
    print(len(real_census))
    print(random_census.shape)
    avg_random_census = []
    z_score = []
    for p in range(len(real_census)):
        N_real_p = real_census[p]
        N_rand_p = np.mean(random_census[:,p])
        std = np.std(random_census[:,p])

        z_p =  ((N_real_p - N_rand_p)/std if std != 0 else 0)
        z_score.append(z_p)
        avg_random_census.append(N_rand_p)

    sp = []
    for i in range(len(z_score)):
        z_norm = np.linalg.norm(z_score)
        norm_z_score = (z_score[i]/z_norm if z_norm != 0 else z_score[i])
        sp.append(round(norm_z_score,4))

    plot_sp(sp, real_census, avg_random_census, over_rep)
    return sp

def plot_sp(sp, real_frq, random_frq, over_rep):
    print("[3.1] >> Motif and anti-motif")
    D = 0.1
    p = 0.01
    thresholds = [r*D for r in random_frq]
    thr_anti   = [r*D*-1 for r in random_frq]
    norm_motif = []
    norm_anti  = []
    for i in range(len(thresholds)):
        norm_m = np.linalg.norm(thresholds)
        norm_a = np.linalg.norm(thr_anti)
        norm_m_score = (thresholds[i]/norm_m if norm_m != 0 else thresholds[i])
        norm_a_score = (thr_anti[i]/norm_a if norm_a != 0 else thr_anti[i])
        norm_motif.append(round(norm_m_score,4))
        norm_anti.append(round(norm_a_score,4))

    print("ov shape:",np.array(over_rep).shape)


    print("[3.1] >> Motif frequencies")
    width = 0.4
    plt.bar(np.arange(len(real_frq)) - (width/2), real_frq, width=width, color="yellowgreen", label="real network")
    plt.bar(np.arange(len(random_frq)) + (width/2), random_frq, width=width, color="coral", label="avg random networks")
    plt.xticks(np.arange(len(real_frq)), labels=[str(i) for i in np.arange(1,len(real_frq)+1)])
    plt.xlabel("motif ID")
    plt.ylabel("frequency")
    plt.grid(axis="y")
    plt.legend()
    plt.title("Class-3 motif frequency")
    plt.show()

    print("[3.1] >> Significance Profile:",sp)
    patterns = [str(i) for i in np.arange(1,14)]
    plt.plot(np.arange(len(sp)), sp, 'o-', label="real net z-scores")
    plt.plot(np.arange(len(norm_motif)), norm_motif, 'd--', color="coral", label="motif")
    plt.plot(np.arange(len(norm_anti)), norm_anti,   'd--', color="red", label="anti-motif")
    plt.yticks(np.arange(-0.7, 0.8, 0.1))
    plt.xticks(np.arange(len(patterns)), labels=patterns)
    plt.ylabel("normalized Z-score")
    plt.xlabel("motif ID")

    plt.title("Network significance profile")
    plt.grid(True)
    plt.legend()
    plt.show()
    return


def p3_1():
    M = compute_adjacency(load_matrix())
    m_3, M_3 = motif3struct_bin(M)

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

    sp = significanceProfile(M, nrand=100)
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
    m_4, M_4 = motif4struct_bin(M)
    print("[3.4] >> Displaying frequencies of 4-node motifs.")
    plt.bar(np.arange(0,200), m_4)
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



if __name__ == '__main__':
    import time

    start = time.time()
    G = load_conn_graph(conn="pdc", freq=10, run="R01")
    p3_1()
    #p3_4()

    end = time.time()
    print("Elapsed time:", (end - start)/60, "min")
