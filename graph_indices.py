from commons import *
from connectivity_graph import load_matrix, compute_adjacency, load_channel_coordinates, load_conn_graph, p1_5


def getKey(item):
    return item[1]


def graph_indices_part_2_1(adj_mat, plots=True, verbose=True):
    # create graph using adjacency matrix
    G_Real = nx.from_numpy_matrix(adj_mat, create_using=nx.DiGraph)
    if verbose:
        print("Resulting network density:", np.sum(adj_mat)/4032)
    
    """
        using density lower than 2% raise exception with using
        average_shortest_path_length in networkx, so we computed
        values for each node and applied averageing based on 
        formulas mebntioned in slides.
    """
    # Global 
    # Cf_real = nx.average_clustering(G_Real)
    # PL_real = nx.average_shortest_path_length(G_Real)
    # print("Resulting global Graph Clustering Coefficient:", Cf_real)
    # print("Resulting global Graph Path Lenght:", PL_real)
    
    # using local data to compute values
    nodes_idx = [i for i in range(64)]  # [0,1,2..., N]
    Cfs_real = list(nx.clustering(G_Real, nodes=nodes_idx).values())
    Cf_avg_real = np.sum(Cfs_real) / 64
    if verbose:
        print("Resulting global Graph Clustering Coefficient using local clustering Coefficient:", Cf_avg_real)
    
    
    PLs_real = nx.shortest_path_length(G_Real)
    sum_path_lenghs = 0
    for ps in PLs_real:
        sum_path_lenghs = sum_path_lenghs + np.sum(list(ps[1].values()))
        
    PL_avg_real = sum_path_lenghs/(4032)
    if verbose:
        print("Resulting global Graph Path Lenght using local shortest path lenghts:", PL_avg_real)
    
    
    degrees     = list(G_Real.degree(nodes_idx))     # degree of all nodes as tuple of (node,degree) form
    in_degrees  = list(G_Real.in_degree(nodes_idx))  # degree of all nodes as tuple of (node,degree) form
    out_degrees = list(G_Real.out_degree(nodes_idx)) # degree of all nodes as tuple of (node,degree) form
    
    
    degree_sorted = sorted(degrees, key=getKey, reverse=False)
    top_10_degrees = degree_sorted[-10:]
    if verbose:
        print("Resulting highest 10 degrees (node, degree) format:", top_10_degrees)

    if plots:
        my_colors = [
            ["blue","red","navy","yellow","cyan","gray","brown","magenta","orange","lime"],
            ["red","green","blue","yellow","cyan","gray","brown","magenta","orange","lime"],
            ["red","green","blue","yellow","cyan","gray","brown","magenta","orange","lime"]
        ]
        plot2_1(degrees,in_degrees,out_degrees,my_colors)

    return Cf_avg_real, PL_avg_real


def plot2_1(node_degs, in_degs, out_degs, colors):
    """
    Plots top 10 nodes for local indices: degree, in-degree, out-degree 
    """
    fig, axs = plt.subplots(1, 3, figsize=(9, 3))
    ch_names = load_channel_coordinates(label=False,map_ch=True)

    degree_sorted      = sorted(node_degs, key=getKey, reverse=False)
    top_10_degrees     = degree_sorted[-10:]
    nodeid_best_degree = [ch_names[n[0]] for n in top_10_degrees]
    best_degree        = [n[1] for n in top_10_degrees]

    axs[0].barh(np.arange(0,10), best_degree, color=colors[0])
    axs[0].set_yticks(np.arange(0,10))
    axs[0].set_yticklabels(nodeid_best_degree)
    axs[0].set_ylabel("channel ID")
    axs[0].set_xlabel("node degree")
    
    in_degree_sorted      = sorted(in_degs, key=getKey, reverse=False)
    top_10_in_degrees     = in_degree_sorted[-10:]
    nodeid_best_in_degree = [ch_names[n[0]] for n in top_10_in_degrees]
    best_in_degree        = [n[1] for n in top_10_in_degrees]
  
    axs[1].barh(np.arange(0,10), best_in_degree, color=colors[1])
    axs[1].set_yticks(np.arange(0,10)),
    axs[1].set_yticklabels(nodeid_best_in_degree)
    axs[1].set_xlabel("in degree")
    
    out_degree_sorted       = sorted(out_degs, key=getKey, reverse=False)
    top_10_out_degrees      = in_degree_sorted[-10:]
    nodeid_best_out_degree = [ch_names[n[0]] for n in top_10_out_degrees]
    best_out_degree         = [n[1] for n in top_10_out_degrees]

    axs[2].barh(np.arange(0,10), best_out_degree, color=colors[2])
    axs[2].set_yticks(np.arange(0,10)),
    axs[2].set_yticklabels(nodeid_best_out_degree)
    axs[2].set_xlabel("out degree")
    
    fig.suptitle("Top 10 channels for local     ")
    plt.show()
    return


def graph_indices_part_2_2(Cf_real, PL_real, random_graph='erdos'):
    if random_graph == 'erdos':
        G_Rand = nx.erdos_renyi_graph(n=64, p=0.4, directed=True)
    else:
        G_Rand = nx.watts_strogatz_graph(n=64, p=0.4, k=5)
        
    Cf_rand = nx.average_clustering(G_Rand)
    PL_rand = nx.average_shortest_path_length(G_Rand)
    
    
    Small_worldness = (Cf_real/Cf_rand)/(PL_real/PL_rand)
    print("Resulting Small worldness:", Small_worldness)
    
    return Small_worldness


def graph_indices_part_2_4(conn_mat, thresholds):
    cl_coeffs = []
    avg_pl    = []
    for threshold in thresholds:
        adj_mat = compute_adjacency(conn_mat, threshold=threshold)  # 0.04597 for PDC
        print("Threshold: ", threshold," Resulting network density:", np.sum(adj_mat)/max_edges)
        
        cl, pl = graph_indices_part_2_1(adj_mat,plots=False)
        print('\n')        

        cl_coeffs.append(cl)
        avg_pl.append(pl)
    
    plot2_4(cl_coeffs,avg_pl)
    return


def plot2_4(cl_coeffs, avg_pl, densities=['1%', '5%', '10%', '20%', '30%', '50%']):
    fig, ax = plt.subplots(1, 2, figsize=(9, 3))
    width = 0.4

    ax[0].bar(np.arange(len(cl_coeffs)), cl_coeffs, width=width, color="yellowgreen",label="avg clustering coefficient")
    ax[0].set_xticks(np.arange(len(densities)))
    ax[0].set_xticklabels(densities)
    ax[0].set_yticks(np.arange(0.0,1.3,0.1))
    ax[0].set_xlabel("network density")
    ax[0].set_ylabel("avg clustering coefficient")
    ax[0].grid(axis="y")
    ax[0].legend()

    ax[1].bar(np.arange(len(avg_pl)), avg_pl, width=0.4, color="coral", label="avg path length")
    ax[1].set_xticks(np.arange(len(densities)))
    ax[1].set_xticklabels(densities)
    ax[1].set_yticks(np.arange(0.0,1.3,0.1))
    ax[1].set_xlabel("network density")
    ax[1].set_ylabel("avg path length")
    ax[1].grid(axis="y")
    ax[1].legend()

    fig.suptitle("Global graph indices per network density")
    plt.show()
    return


def graph_indices_part_2_7(conn_mat, threshold):
    weights = np.random.uniform(0, 2, (64,64))  # random generated weght with shape (n,n)
    
    weights = conn_mat
    
    adj_mat = compute_adjacency(conn_mat, threshold=threshold)  # 0.04597 for PDC
    print("Resulting network density:", np.sum(adj_mat)/max_edges)
    
    
    weighted_adj_mat = np.multiply(weights, adj_mat) # mask weiths using adjacency matrix
    
    G_weighted = nx.from_numpy_matrix(weighted_adj_mat, create_using=nx.DiGraph)
    #print(G_weighted.edges(data="weight"))

    for u, v, weight in G_weighted.edges(data="weight"):
        if weight is not None:
            # print(f'weight ==== {u} --> {v}: {weight}')
            pass
        else:
            print('no weight ==== {u} --> {v} ')
            

    Cfs_wght = list(nx.clustering(G_weighted, weight="weight").values())
    Cf_avg_wght = np.sum(Cfs_wght) / 64
    print("Resulting global Graph Clustering Coefficient using local clustering Coefficient:", Cf_avg_wght)
    
    
    PLs_wght = nx.shortest_path_length(G_weighted, weight="weight")
    sum_path_lenghs = 0
    for ps in PLs_wght:
        # print(list(ps[1].values()))
        sum_path_lenghs += np.sum(list(ps[1].values()))
        
    PL_avg_wght = sum_path_lenghs/(4032)
    print("Resulting global Graph Path Lenght using local shortest path lenghts:", PL_avg_wght)
    
    
    # degrees
    degrees_weighted = list(G_weighted.degree(nodes_idx))  # degree of all nodes as tuple of (node,degree) form
    in_degrees_weighted = list(G_weighted.in_degree(nodes_idx))  # degree of all nodes as tuple of (node,degree) form
    out_degrees_weighted = list(G_weighted.out_degree(nodes_idx))  # degree of all nodes as tuple of (node,degree) form
    
    
    degree_sorted_weighted = sorted(degrees_weighted, key=getKey)
    
    top_10_degrees_weighted = degree_sorted_weighted[-10:]
    print("Resulting highest 10 degrees (node, degree) format:", top_10_degrees_weighted) 

    my_colors = [
        ["blue","red","navy","yellow","cyan","gray","brown","magenta","orange","lime"],
        ["red","green","blue","yellow","cyan","gray","brown","magenta","orange","lime"],
        ["red","green","blue","yellow","cyan","gray","brown","magenta","orange","lime"]
    ]
    plot2_1(degrees_weighted,in_degrees_weighted,out_degrees_weighted,colors=my_colors)

    return  


def p2_3(freq, run, threshold_pdc, threshold_dtf):
    """
    Computes for both pdc and dtf method:
        - cl (Average Clustering Coefficient)
        - pl (Average Path Length)
    """
    pdc_mat = compute_adjacency(load_matrix(conn_method='pdc', freq=freq, run=run), threshold=threshold_pdc)    
    dtf_mat = compute_adjacency(load_matrix(conn_method='dtf', freq=freq, run=run), threshold=threshold_dtf)
    cl_pdc, pl_pdc = graph_indices_part_2_1(pdc_mat, plots=False, verbose=False)
    cl_dtf, pl_dtf = graph_indices_part_2_1(dtf_mat, plots=False, verbose=False)

    print("[2.3] >> Average Clustering Coefficient PDC: {:.2f}%".format(100*cl_pdc))
    print("[2.3] >> Average Clustering Coefficient DTF: {:.2f}%".format(100*cl_dtf))
    print("[2.3] >> Average PAth Length PDC: {:.2f}".format(pl_pdc))
    print("[2.3] >> Average PAth Length DTF: {:.2f}".format(pl_dtf))


def p2_5(G):
    print("[2.5] >>")
    p1_5(G, point='2.5')


def p2_6(run):
    """
    Compares the graph indices found in p2_1 with those of another band's frequeny\n
    Thresold values are stored here and not calculated again
    """
    if run == 'R01':
        threshold_10Hz = THRES_10HZ_R01_20percent
        threshold_25Hz = 0.1268000000000301
    elif run == 'R02':
        threshold_10Hz = THRES_10HZ_R02_20percent
        threshold_25Hz = 0.12200000000003064
    
    pdc_mat_10Hz = compute_adjacency(load_matrix(conn_method='pdc', freq=10, run=run), threshold=threshold_10Hz)    
    cl_pdc_10Hz, pl_pdc_10Hz = graph_indices_part_2_1(pdc_mat_10Hz, plots=False, verbose=False)

    pdc_mat_25Hz = compute_adjacency(load_matrix(conn_method='pdc', freq=25, run=run), threshold=threshold_25Hz)    
    cl_pdc_25Hz, pl_pdc_25Hz = graph_indices_part_2_1(pdc_mat_25Hz,  plots=False, verbose=False)

    print("[2.3] >> Average Clustering Coefficient PDC - 10Hz: {:.4f}%".format(100*cl_pdc_10Hz))
    print("[2.3] >> Average Clustering Coefficient PDC - 25Hz: {:.4f}%".format(100*cl_pdc_25Hz))
    print("[2.3] >> Average PAth Length PDC - 10Hz: {:.4f}".format(pl_pdc_10Hz))
    print("[2.3] >> Average PAth Length PDC - 25Hz: {:.4f}".format(pl_pdc_25Hz))


def p2_1(conn, freq, run):
    if run == 'R01':
        threshold = THRES_10HZ_R01_20percent
    elif run == 'R02':
        threshold = THRES_10HZ_R02_20percent

    adj_mat = compute_adjacency(load_matrix(conn_method=conn, freq=freq, run=run), threshold=threshold)
    Cf_real, PL_real = graph_indices_part_2_1(adj_mat)



if __name__ == '__main__':
    conn_method = 'dtf'
    random_graph = 'erdos' # 'watts' and erdos are options
    N = 64
    
    max_edges = 4032
    '''
    # get adjacency matrix
    conn_mat = load_matrix(conn_method=conn_method)
    print("mat shape:", conn_mat.shape)

    threshold_20_percent_density = 0.137 # 0.04597 for PDC
    adj_mat = compute_adjacency(conn_mat, threshold=threshold_20_percent_density)  
    print("Resulting network density:", np.sum(adj_mat)/max_edges)

    print('\n================== P 2.1 ==============================')
    Cf_real, PL_real = graph_indices_part_2_1(adj_mat)


    print('\n================== P 2.2 ==============================')
    print('small worls formula = (Cf_G/Cf_rand)/(PL_G/PL_rand)')
    # small worls formula = (Cf_G/Cf_rand)/(PL_G/PL_rand)
    Small_worldness = graph_indices_part_2_2(Cf_real, PL_real, random_graph)



    print('\n================== P 2.4 ==============================')
    # just change threshold values in crearting adjacency matrix to tune density as 
    # mentioned in P 1.3
    #  densities = [1%, 5%, 10%, 20%, 30%, 50%]
    thresholds = [0.41, 0.24, 0.187, 0.137, 0.1, 0.055]
    graph_indices_part_2_4(conn_mat, thresholds)



    print('\n================== P 2.7 ==============================')
    threshold = threshold_20_percent_density
    graph_indices_part_2_7(conn_mat, threshold)
    '''


    print('\n================== P 2.3 ==============================')
    p2_3(freq=10, run='R01', threshold_pdc=0.1226, threshold_dtf=0.1378)

    #print('\n================== P 2.5 ==============================')
    #G = load_conn_graph(conn="pdc", freq=10, run="R01")
    #p2_5(G)

    print('\n================== P 2.6 ==============================')
    p2_6('R01')
    p2_6('R02')
