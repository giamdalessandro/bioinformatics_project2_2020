import louvain as lv
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import networkx as nx
import igraph as ig

from connectivity_graph import load_conn_graph, load_channel_coordinates, load_matrix, compute_adjacency

conn_mat = load_matrix(conn_method="pdc", freq=10, run="R01", auto='auto')
adj_mat  = compute_adjacency(conn_mat, threshold=0.1226)

G = ig.Graph.Adjacency((adj_mat > 0).tolist())

ig.summary(G)

with open("data/channel_locations.txt") as f:
    mapping = []
    for line in f:
        l = line.split(sep='        ')
        if l[0] != '\ufeff#':
            mapping.append(str(l[1]))

G.vs['label'] = mapping

ig.summary(G)








'''
    G = load_conn_graph(conn='pdc', freq=10, run='R01')
    print("Graph has {} nodes and {} edges".format(len(G.nodes()), len(G.edges())))

    pos = load_channel_coordinates()

    dendo = cl.generate_dendrogram(G)
    print("len dendogram", len(dendo))

    for level in range(len(dendo)-1) :
        partition = cl.partition_at_level(dendo, level)
        print("modularity = ", cl.modularity(partition, G))

        # color the nodes according to their partition
        cmap = cm.get_cmap('viridis', max(partition.values()) + 1)
        nx.draw_networkx_nodes(G, pos, partition.keys(), node_size=40, cmap=cmap, node_color=list(partition.values()))
        nx.draw_networkx_edges(G, pos, alpha=0.5)
        plt.show()

    # compute the best partition
    partition = cl.best_partition(G)
    # color the nodes according to their partition
    cmap = cm.get_cmap('viridis', max(partition.values()) + 1)
    nx.draw_networkx_nodes(G, pos, partition.keys(
    ), node_size=40, cmap=cmap, node_color=list(partition.values()))
    nx.draw_networkx_edges(G, pos, alpha=0.5)
    plt.show()
    print("best modularity = ", cl.modularity(partition, G))
'''
