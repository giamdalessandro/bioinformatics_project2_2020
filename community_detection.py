import matplotlib.cm as cm
import matplotlib.pyplot as plt
import networkx as nx
import community as cl

from connectivity_graph import load_conn_graph, load_channel_coordinates, load_matrix, compute_adjacency, p1_5


def relabel_partition(partition):
    """
    Relabels node in the partition obtained using igraph\n
    It's not optimized nor elegant, but the graph is small, so...
    """
    print('\n', partition, '\n')
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
        
    partition = {}
    print("[4.1] >> Partitions found: {}\n".format(len(d)))
    for k in d.keys():
        print("[Partition {}]".format(k))
        for p in d[k]:
            partition.update({ p : k})
            print(p, end=" ", flush=True)
        print('\n')
    return partition


# NOTE: to check, dunno if it's really the best partition
def best_partition_louvain(conn_method="pdc", freq=10, run="R01", auto='auto', threshold=0.1226):
    """
    Since there's no library for computing the best partition using the Louvain algorithm on directed graphs for\n
    `NetworkX`, we rely on `iGraph` only for this function.
    """
    import louvain
    import igraph as ig
    conn_mat = load_matrix(conn_method=conn_method, freq=freq, run=run, auto=auto)
    adj_mat  = compute_adjacency(conn_mat, threshold=threshold)
    G = ig.Graph.Adjacency((adj_mat > 0).tolist())
    partition = louvain.find_partition(G, louvain.ModularityVertexPartition)   # finds best partition
    return partition
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


def p4_1(conn_method="pdc", freq=10, run="R01", auto='auto', threshold=0.1226):
    
    partition = best_partition_louvain(conn_method=conn_method, freq=freq,
                                       run=run, auto=auto, threshold=threshold)
    return relabel_partition(partition)


def p4_2(G, partition, algorithm='Louvain'):
    """
    Display a topographical representation of the community structure
    """
    pos = load_channel_coordinates()
    cmap = cm.get_cmap('viridis', max(partition.values()) + 1)
    nx.draw_networkx(G, pos=pos, arrows=True, with_labels=True, nodelist=partition.keys(), node_size=700,
                     cmap=cmap, node_color=list(partition.values()), edge_color='black')
    plt.title("Topographical representation of the community structure found using the {} algorithm".format(algorithm))
    plt.show()


G = load_conn_graph(conn="pdc", freq=10, run="R01",
                    auto='auto', threshold=0.1226)
print("Graph has {} nodes and {} edges".format(len(G.nodes()), len(G.edges())))

partition = p4_1()
p4_2(G, partition)


