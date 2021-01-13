import community as cl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import networkx as nx

# load the karate club graph
G = nx.random_lobster(100, 0.9, 0.9)
print("Graph has {} nodes and {} edges".format(len(G.nodes()), len(G.edges())))

dendo = cl.generate_dendrogram(G)
print("len dendogram", len(dendo))
for level in range(len(dendo)-1) :
    partition = cl.partition_at_level(dendo, level)
    pos = nx.circular_layout(G)
    # color the nodes according to their partition
    cmap = cm.get_cmap('viridis', max(partition.values()) + 1)
    nx.draw_networkx_nodes(G, pos, partition.keys(), node_size=40, cmap=cmap, node_color=list(partition.values()))
    nx.draw_networkx_edges(G, pos, alpha=0.5)
    plt.show()
    print("modularity = ", cl.modularity(partition, G))

# compute the best partition
partition = cl.best_partition(G)
pos = nx.circular_layout(G)
# color the nodes according to their partition
cmap = cm.get_cmap('viridis', max(partition.values()) + 1)
nx.draw_networkx_nodes(G, pos, partition.keys(
), node_size=40, cmap=cmap, node_color=list(partition.values()))
nx.draw_networkx_edges(G, pos, alpha=0.5)
plt.show()
print("best modularity = ", cl.modularity(partition, G))
