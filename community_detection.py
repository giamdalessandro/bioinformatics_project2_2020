import matplotlib.cm as cm
import matplotlib.pyplot as plt
import networkx as nx
import community as cl
import infomap

from commons import *
from connectivity_graph import load_conn_graph, load_channel_coordinates, load_matrix, compute_adjacency, p1_5, map_index_to_channels


def relabel_partition_louvain(partition):
    """
    Relabels node in the partition obtained using igraph\n
    It's not optimized nor elegant, but the graph is small, so...
    """
    partition = list(partition)
    d = {}

    mapping = map_index_to_channels()
    for i in range(len(partition)):
        part = []
        for p in partition[i]:
            part.append(mapping[p])
        d.update({i: part})
    partition = {}
    print("[4.1] >> Communities found using Louvain algorithm: {}\n".format(len(d)))
    for k in d.keys():
        print("[Community {}]".format(k))
        for p in d[k]:
            partition.update({ p : k})
            print(p, end=" ", flush=True)
        print('\n')
    return partition, d


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
    

def best_partition_infomap(G):
    """
    Partition network with the Infomap algorithm.
    Annotates nodes with 'community' id.
    """

    im = infomap.Infomap("--directed")  # --prefer-modular-solution yields more cluster but with less accuracy

    print("Building Infomap network from a NetworkX graph...")
    for source, target in G.edges:
        im.add_link(source, target)

    print("Find communities with Infomap...")
    im.run()

    print(f"Found {im.num_top_modules} modules with codelength: {im.codelength}")

    communities = im.get_modules()
    nx.set_node_attributes(G, communities, 'community')
    return communities


####


def p4_1(conn_method="pdc", freq=10, run="R01", auto='auto'):
    if run == 'R01':
        threshold = THRES_PDC_10HZ_R01_20percent
    elif run == 'R02':
        threshold = THRES_PDC_10HZ_R02_20percent

    partition = best_partition_louvain(conn_method=conn_method, freq=freq, run=run, auto=auto, threshold=threshold)
    return relabel_partition_louvain(partition)


def p4_2(run, partition):
    """
    Display a topographical representation of the community structure
    """
    if run == 'R01':
        threshold = THRES_PDC_10HZ_R01_20percent
    elif run == 'R02':
        threshold = THRES_PDC_10HZ_R02_20percent
    G = load_conn_graph(conn="pdc", freq=10, run=run, auto="auto", threshold=threshold)
    p1_5(G, point='4.2', communities=partition)


def p4_3(conn_method="pdc", freq=10, run='R01', auto='auto'):
    if run == 'R01':
        threshold = THRES_PDC_10HZ_R01_20percent
    elif run == 'R02':
        threshold = THRES_PDC_10HZ_R02_20percent
    conn_mat = load_matrix(conn_method=conn_method, freq=freq, run=run, auto=auto)
    adj_mat = compute_adjacency(conn_mat, threshold=threshold)
    G = nx.from_numpy_array(adj_mat, create_using=nx.DiGraph)

    communities = best_partition_infomap(G)
    mapping = map_index_to_channels()
    G = nx.relabel_nodes(G, mapping)

    ### this is just to print the result, not needed for the logic
    l = []
    for i in range(len(set(communities.values()))):
        l.append(list())
    for node, comm in communities.items():
        l[comm-1].append(mapping[node])
    print("[4.3] >> Communities found using Infomap algorithm: {}\n".format(len(l)))
    for i in range(len(l)):
        print("[Community {}]".format(i))
        for node in l[i]:
            print(node, end=" ", flush=True)
        print('\n')
    ###

    communities = [c - 1 for c in nx.get_node_attributes(G, 'community').values()]
    p1_5(G, point='4.3', communities=communities)
    return l


def jaccard(S1, S2):
    S1 = set(S1)
    S2 = set(S2)
    return 100*len(S1.intersection(S2))/len(S1.union(S2))



### main

if __name__ == '__main__':
    G = load_conn_graph(conn="pdc", freq=10, run="R01", auto='auto', threshold=0.1226)
    print("Graph has {} nodes and {} edges".format(len(G.nodes()), len(G.edges())))

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
