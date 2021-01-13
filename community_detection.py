import louvain
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import networkx as nx
import igraph as ig
import community as cl

from connectivity_graph import load_conn_graph, load_channel_coordinates, load_matrix, compute_adjacency, p1_5


def relabel_partition(partition):
    """
    Relabels node in the partition obtained using igraph
    """
    partition = list(partition)
    d = {}

    def map_index_to_channels():
        """
        Maps channel coordinates to indeces in the file
        """
        with open("data/channel_locations.txt") as f:
            pos = {}
            for line in f:
                # yes, there are 8 spaces in the file.
                l = line.split(sep='        ')
                if l[0] != '\ufeff#':
                    pos.update({int(l[0])-1: str(l[1])})
        return pos
    
    mapping = map_index_to_channels()
    for i in range(len(partition)):
        part = []
        for p in partition[i]:
            part.append(mapping[p])
        d.update({i: part})
    return d


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

partition, diff = louvain.find_partition(G, louvain.ModularityVertexPartition)   # finds best partition
'''
    best_partition = louvain.find_partition(G, louvain.CPMVertexPartition, resolution_parameter=0.1)
    print(best_partition)
    #ig.plot(best_partition, layout=G.layout("kk"))


    partition = louvain.CPMVertexPartition(G, resolution_parameter=1.0)
    optimiser = louvain.Optimiser()
    old = optimiser.optimise_partition(partition)
    print("len partition", len(partition), '\tdiff =', old)

    rc = 0.9
    while True:
        partition = louvain.CPMVertexPartition(G, resolution_parameter=rc)
        optimiser = louvain.Optimiser()
        diff = optimiser.optimise_partition(partition)
        print("len partition", len(partition), '\tdiff =', diff-old)
        if diff - old < 0 or len(partition) == 1:
            break
        old = diff
        rc -= 0.1
'''



partition = relabel_partition(partition)

G = load_conn_graph(conn='pdc', freq=10, run='R01')
print("Graph has {} nodes and {} edges".format(len(G.nodes()), len(G.edges())))


pos = load_channel_coordinates()
cmap = cm.get_cmap('viridis', max(partition.values()) + 1)
nx.draw_networkx_nodes(G, pos=pos, nodelist=partition.keys(), node_size=40,
                       cmap=cmap, node_color=list(partition.values()))
nx.draw_networkx_edges(G, pos, alpha=0.5)
plt.show()


# G = nx.erdos_renyi_graph(100, 0.01)
# pos = nx.spring_layout(G)
# partition = cl.best_partition(G)
# 
# print(partition)
# 
# cmap = cm.get_cmap('viridis', max(partition.values()) + 1)
# 
# nx.draw_networkx_nodes(G, pos, partition.keys(), node_size=40, cmap=cmap, node_color=list(partition.values()))
# 
# nx.draw_networkx_edges(G, pos, alpha=0.5)
# 
# plt.show()
# print("best modularity = ", cl.modularity(partition, G))

