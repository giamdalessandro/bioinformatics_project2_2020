import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pathlib

import infomap

from connectivity_graph import load_matrix, compute_adjacency, p1_5


def find_communities(G):
    """
    Partition network with the Infomap algorithm.
    Annotates nodes with 'community' id.
    """

    im = infomap.Infomap("--directed --prefer-modular-solution")

    print("Building Infomap network from a NetworkX graph...")
    for source, target in G.edges:
        im.add_link(source, target)

    print("Find communities with Infomap...")
    im.run()

    print(f"Found {im.num_top_modules} modules with codelength: {im.codelength}")

    communities = im.get_modules()
    print("[4.3] >> Communities found;", communities)
    nx.set_node_attributes(G, communities, 'community')


def draw_network(G):
    # position map
    pos = nx.spring_layout(G)
    # community index
    communities = [c - 1 for c in nx.get_node_attributes(G, 'community').values()]
    num_communities = max(communities) + 1

    # color map from http://colorbrewer2.org/
    cmap_light = colors.ListedColormap(
        ['#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f', '#cab2d6'], 'indexed', num_communities)
    cmap_dark = colors.ListedColormap(
        ['#1f78b4', '#33a02c', '#e31a1c', '#ff7f00', '#6a3d9a'], 'indexed', num_communities)

    # edges
    nx.draw_networkx_edges(G, pos)

    # nodes
    node_collection = nx.draw_networkx_nodes(
        G, pos=pos, node_color=communities, cmap=cmap_light)

    # set node border color to the darker shade
    dark_colors = [cmap_dark(v) for v in communities]
    node_collection.set_edgecolor(dark_colors)

    # Print node labels separately instead
    for n in G.nodes:
        plt.annotate(n,
                     xy=pos[n],
                     textcoords='offset points',
                     horizontalalignment='center',
                     verticalalignment='center',
                     xytext=[0, 2],
                     color=cmap_dark(communities[n]))

    plt.axis('off')
    plt.show()


conn_mat = load_matrix(conn_method="pdc", freq=10, run='R01', auto='auto')
adj_mat = compute_adjacency(conn_mat, threshold=0.1226)
G = nx.from_numpy_array(adj_mat, create_using=nx.DiGraph)
print("Graph has {} nodes and {} edges".format(len(G.nodes()), len(G.edges())))

find_communities(G)
with open("data/channel_locations.txt") as f:
    mapping = {}
    for line in f:
        l = line.split(sep='        ')
        if l[0] != '\ufeff#':
            mapping.update({int(l[0]) - 1: str(l[1])})
G = nx.relabel_nodes(G, mapping)
#draw_network(G)
communities = [c - 1 for c in nx.get_node_attributes(G, 'community').values()]
p1_5(G, point='4.3', communities=communities)
