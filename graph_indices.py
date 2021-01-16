from commons import *
from connectivity_graph import load_matrix, compute_adjacency, load_channel_coordinates, load_conn_graph, p1_5

from random import randint

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
        my_colors = plt.get_cmap('tab20')
        # my_colors = [
        #     ["blue","red","navy","yellow","cyan","gray","brown","magenta","orange","lime"],
        #     ["red","green","blue","yellow","cyan","gray","brown","magenta","orange","lime"],
        #     ["red","green","blue","yellow","cyan","gray","brown","magenta","orange","lime"]
        # ]
        plot2_1(degrees,in_degrees,out_degrees,my_colors)

    return Cf_avg_real, PL_avg_real

def plot2_1(node_degs, in_degs, out_degs, colors):
    """
    Plots top 10 nodes for local indices: degree, in-degree, out-degree 
    """
    rescale = lambda y: (y - np.min(y)) / (np.max(y) - np.min(y))

    #fig, axs = plt.subplots(1, 2, figsize=(9, 3))
    ch_names = load_channel_coordinates(label=False,map_ch=True)

    degree_sorted      = sorted(node_degs, key=getKey, reverse=False)
    top_10_degrees     = degree_sorted[-10:]
    nodeid_best_degree = [ch_names[n[0]] for n in top_10_degrees]
    best_degree        = [n[1] for n in top_10_degrees]

    #axs[0].barh(np.arange(0, 10), best_degree, color=colors(rescale(best_degree)))
    #axs[0].set_yticks(np.arange(0,10))
    #axs[0].set_yticklabels(nodeid_best_degree)
    #axs[0].set_ylabel("channel ID")
    #axs[0].set_xlabel("node degree")
    #axs[0].grid()
    
    in_degree_sorted      = sorted(in_degs, key=getKey, reverse=False)
    top_10_in_degrees     = in_degree_sorted[-10:]
    nodeid_best_in_degree = [ch_names[n[0]] for n in top_10_in_degrees]
    best_in_degree        = [n[1] for n in top_10_in_degrees]
  
    #axs[1].barh(np.arange(0, 10), best_in_degree, color=colors(rescale(best_in_degree)))
    #axs[1].set_yticks(np.arange(0,10)),
    #axs[1].set_yticklabels(nodeid_best_in_degree)
    #axs[1].set_xlabel("in degree")
    #axs[1].grid()
    #
    out_degree_sorted       = sorted(out_degs, key=getKey, reverse=False)
    top_10_out_degrees      = out_degree_sorted[-10:]
    nodeid_best_out_degree = [ch_names[n[0]] for n in top_10_out_degrees]
    best_out_degree         = [n[1] for n in top_10_out_degrees]

    #axs[2].barh(np.arange(0,10), best_out_degree, color=colors(rescale(best_out_degree)))
    #axs[2].set_yticks(np.arange(0,10)),
    #axs[2].set_yticklabels(nodeid_best_out_degree)
    #axs[2].set_xlabel("out degree")
    #axs[2].grid()
    #

    plt.scatter(nodeid_best_degree, best_degree,        s=100, label="10 best degree nodes", cmap=colors)
    plt.scatter(nodeid_best_in_degree, best_in_degree,  s=100, label="10 best in degree nodes", cmap=colors)
    plt.scatter(nodeid_best_out_degree, best_out_degree,s=100, label="10 best out degree nodes", cmap=colors)
    plt.grid()
    plt.legend()
    #fig.suptitle("Top 10 channels for local     ")
    plt.show()
    return


def graph_indices_part_2_2(cf_real, pl_real, iters=500, cf_rand=None, pl_rand=None):
    if cf_rand is not None and pl_rand is not None:
        return (cf_real/cf_rand)/(pl_real/pl_rand)
    cf_rand = 0
    pl_rand = 0
    for i in range(iters):
        #G_Rand = nx.erdos_renyi_graph(n=64, p=randint(4,6)/10, seed=i, directed=True)        
        G_Rand = nx.erdos_renyi_graph(n=64, p=0.4, seed=i, directed=True)        
        cf = nx.average_clustering(G_Rand)
        pl = nx.average_shortest_path_length(G_Rand)
        cf_rand += cf 
        pl_rand += pl
    
    cf_rand = cf_rand / iters
    pl_rand = pl_rand / iters
    print("cf real / cf rand = ", cf_real, '/', cf_rand)
    print("pl real / l rand = ", pl_real, '/', pl_rand)
    Small_worldness = (cf_real/cf_rand)/(pl_real/pl_rand)    
    return Small_worldness


def cf_pl(iters=500):
    cf_rand = 0
    pl_rand = 0
    for i in range(iters):
        #G_Rand = nx.erdos_renyi_graph(n=64, p=randint(4,6)/10, seed=i, directed=True)
        G_Rand = nx.erdos_renyi_graph(n=64, p=0.4, seed=i, directed=True)
        cf = nx.average_clustering(G_Rand)
        pl = nx.average_shortest_path_length(G_Rand)
        cf_rand += cf
        pl_rand += pl

    cf_rand = cf_rand / iters
    pl_rand = pl_rand / iters
    print("cf_rand = ", cf_rand)
    print("pl_rand = ", pl_rand)
    return cf_rand, pl_rand


def graph_indices_part_2_4(conn_mat, thresholds, cf_rand=None, pl_rand=None):
    cl_coeffs = []
    avg_pl    = []
    smalls    = []
    for threshold in thresholds:
        # if threshold == THRES_PDC_10HZ_R01_01percent or threshold == THRES_PDC_10HZ_R02_01percent:
        #     d = 0.1
        # elif threshold == THRES_PDC_10HZ_R01_05percent or threshold == THRES_PDC_10HZ_R01_05percent:
        #     d = 0.1
        # elif threshold == THRES_PDC_10HZ_R01_10percent or threshold == THRES_PDC_10HZ_R01_10percent:
        #     d = 0.1
        # elif threshold == THRES_PDC_10HZ_R01_20percent or threshold == THRES_PDC_10HZ_R01_20percent:
        #     d = 0.2
        # elif threshold == THRES_PDC_10HZ_R01_30percent or threshold == THRES_PDC_10HZ_R01_30percent:
        #     d = 0.3
        # elif threshold == THRES_PDC_10HZ_R01_50percent or threshold == THRES_PDC_10HZ_R01_50percent:
        #     d = 0.5
        adj_mat = compute_adjacency(conn_mat, threshold=threshold)        
        cl, pl = graph_indices_part_2_1(adj_mat, plots=False, verbose=False)
        small_worldness = graph_indices_part_2_2(cl, pl, cf_rand=cf_rand, pl_rand=pl_rand) #, density=d)
        cl_coeffs.append(cl)
        avg_pl.append(pl)
        smalls.append(small_worldness)
    
    plot2_4(cl_coeffs,avg_pl,smalls)
    return


def plot2_4(cl_coeffs, avg_pl, smalls, densities=['1%', '5%', '10%', '20%', '30%', '50%']):
    fig, ax = plt.subplots(1, 3, figsize=(8, 4))
    width = 0.4
    #max_y = max(avg_pl) if max(avg_pl) >= max(cl_coeffs) else max(cl_coeffs)

    max_y = max(max(cl_coeffs), max(avg_pl))

    ax[0].bar(np.arange(len(cl_coeffs)), cl_coeffs, width=width, color="yellowgreen",label="avg clustering coefficient")
    ax[0].set_xticks(np.arange(len(densities)))
    ax[0].set_xticklabels(densities)
    ax[0].set_yticks(np.arange(0.0,max_y,0.2))
    ax[0].set_xlabel("network density")
    ax[0].set_ylabel("avg clustering coefficient")
    ax[0].grid(axis="y")
    ax[0].legend()

    ax[1].bar(np.arange(len(avg_pl)), avg_pl, width=0.4, color="coral", label="avg path length")
    ax[1].set_xticks(np.arange(len(densities)))
    ax[1].set_xticklabels(densities)
    ax[1].set_yticks(np.arange(0.0,max_y,0.2))
    ax[1].set_xlabel("network density")
    ax[1].set_ylabel("avg path length")
    ax[1].grid(axis="y")
    ax[1].legend()

    ax[2].bar(np.arange(len(smalls)), smalls, width=0.4, color="darkturquoise", label="small-worldness")
    ax[2].set_xticks(np.arange(len(densities)))
    ax[2].set_xticklabels(densities)
    ax[2].set_yticks(np.arange(0.0, max(smalls) , 0.5))
    ax[2].set_xlabel("network density")
    ax[2].set_ylabel("small-worldness")
    ax[2].grid(axis="y")
    ax[2].legend()

    fig.suptitle("Global graph indices per network density")
    plt.show()
    return


def graph_indices_part_2_7(conn_mat, threshold):
    weights = np.random.uniform(0, 2, (64,64))  # random generated weght with shape (n,n)
    weights = conn_mat # ???
    adj_mat = compute_adjacency(conn_mat, threshold=threshold)  # 0.04597 for PDC    
    
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
    # print("Resulting global Graph Clustering Coefficient using local clustering Coefficient:", Cf_avg_wght)
    
    
    PLs_wght = nx.shortest_path_length(G_weighted, weight="weight")
    sum_path_lenghs = 0
    for ps in PLs_wght:
        # print(list(ps[1].values()))
        sum_path_lenghs += np.sum(list(ps[1].values()))
        
    PL_avg_wght = sum_path_lenghs/(4032)
    # print("Resulting global Graph Path Lenght using local shortest path lenghts:", PL_avg_wght)
    
    
    # degrees
    nodes_idx = [i for i in range(64)]  # [0,1,2..., N]
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

    return Cf_avg_wght, PL_avg_wght


def p2_1(conn, freq, run):
    if run == 'R01':
        threshold = THRES_PDC_10HZ_R01_20percent
    elif run == 'R02':
        threshold = THRES_PDC_10HZ_R02_20percent
    print("\n[2.1] >> analyzing run", run)
    adj_mat = compute_adjacency(load_matrix(conn_method=conn, freq=freq, run=run), threshold=threshold)
    cf_real, pl_real = graph_indices_part_2_1(adj_mat, verbose=False)
    print("[2.1] >> Average Clustering Coefficient PDC: {:.6f}".format(cf_real))
    print("[2.1] >> Average Path Length PDC: {:.6f}".format(pl_real))
    return cf_real, pl_real


def p2_2(cf_real, pl_real):
    print('\n[2.2] >> small worlds formula = (Cf_G/Cf_rand)/(PL_G/PL_rand)')
    # small worls formula = (Cf_G/Cf_rand)/(PL_G/PL_rand)
    small_worldness = graph_indices_part_2_2(cf_real, pl_real)   
    print("[2.2] >> Small worldness PDC: {:.6f}".format(small_worldness))
    return small_worldness


def p2_3(freq, run):
    """
    Computes p2_1 for dtf method:
        - cl (Average Clustering Coefficient)
        - pl (Average Path Length)
        - SMALL WORLDNESS
    """
    if run == 'R01':
        threshold_pdc = THRES_PDC_10HZ_R01_20percent
        threshold_dtf = THRES_DTF_10HZ_R01_20percent
    elif run == 'R02':
        threshold_pdc = THRES_PDC_10HZ_R02_20percent
        threshold_dtf = THRES_DTF_10HZ_R02_20percent

    dtf_mat = compute_adjacency(load_matrix(conn_method='dtf', freq=freq, run=run), threshold=threshold_dtf)
    cl_dtf, pl_dtf = graph_indices_part_2_1(dtf_mat, plots=False, verbose=False)
    small_worldness = graph_indices_part_2_2(cl_dtf, cl_dtf)
    print("\n")
    print("[2.3] >> Average Clustering Coefficient DTF: {:.6f}".format(cl_dtf))
    print("[2.3] >> Average Path Length DTF: {:.6f}".format(pl_dtf))
    print("[2.3] >> Small worldness DTF:  {:.6f}".format(small_worldness))





def p2_4(run, cf_rand=None, pl_rand=None):
    # just change threshold values in crearting adjacency matrix to tune density as
    # mentioned in P 1.3
    print('\n[2.4] >> Behaviours of global graph indices in function of network density')
    conn_mat = load_matrix('pdc', 10, run, verbose=False)
    if run == 'R01':
        thresholds = [THRES_PDC_10HZ_R01_01percent,
                      THRES_PDC_10HZ_R01_05percent,
                      THRES_PDC_10HZ_R01_10percent,
                      THRES_PDC_10HZ_R01_20percent,
                      THRES_PDC_10HZ_R01_30percent,
                      THRES_PDC_10HZ_R01_50percent]
    else:
        thresholds = [THRES_PDC_10HZ_R02_01percent,
                      THRES_PDC_10HZ_R02_05percent,
                      THRES_PDC_10HZ_R02_10percent,
                      THRES_PDC_10HZ_R02_20percent,
                      THRES_PDC_10HZ_R02_30percent,
                      THRES_PDC_10HZ_R02_50percent]
    graph_indices_part_2_4(conn_mat, thresholds, cf_rand=cf_rand, pl_rand=pl_rand)


def p2_5(G):
    print("[2.5] >> Displaying topological graph...")
    p1_5(G, point='2.5')


def p2_6(run):
    """
    Compares the graph indices found in p2_1 with those of another band's frequeny\n
    Thresold values are stored here and not calculated again
    """
    if run == 'R01':
        threshold_10Hz = THRES_PDC_10HZ_R01_20percent
        threshold_25Hz = THRES_PDC_25HZ_R01_20percent
    elif run == 'R02':
        threshold_10Hz = THRES_PDC_10HZ_R02_20percent
        threshold_25Hz = THRES_PDC_25HZ_R02_20percent
    
    pdc_mat_10Hz = compute_adjacency(load_matrix(conn_method='pdc', freq=10, run=run), threshold=threshold_10Hz)    
    cl_pdc_10Hz, pl_pdc_10Hz = graph_indices_part_2_1(pdc_mat_10Hz, plots=False, verbose=False)

    pdc_mat_25Hz = compute_adjacency(load_matrix(conn_method='pdc', freq=25, run=run), threshold=threshold_25Hz)    
    cl_pdc_25Hz, pl_pdc_25Hz = graph_indices_part_2_1(pdc_mat_25Hz,  plots=False, verbose=False)

    print("[2.6] >> Average Clustering Coefficient PDC - 10Hz: {:.6f}%".format(cl_pdc_10Hz))
    print("[2.6] >> Average Clustering Coefficient PDC - 25Hz: {:.6f}%".format(cl_pdc_25Hz))
    print("[2.6] >> Average Path Length PDC - 10Hz: {:.6f}".format(pl_pdc_10Hz))
    print("[2.6] >> Average Path Length PDC - 25Hz: {:.6f}".format(pl_pdc_25Hz))
    p2_1('pdc', 25, run)


def p2_7(run):
    if run == 'R01':
        threshold = THRES_PDC_10HZ_R01_20percent
    elif run == 'R02':
        threshold = THRES_PDC_10HZ_R02_20percent
    conn_mat = load_matrix('pdc', 10, run, verbose=False)
    cl, pl = graph_indices_part_2_7(conn_mat, threshold)
    print("[2.7] >> Average Weighted Clustering Coefficient PDC: {:.6f}".format(cl))
    print("[2.7] >> Average Weighted Path Length PDC: {:.6f}".format(pl))

if __name__ == '__main__':
    conn_method = 'dtf'
    random_graph = 'erdos' # 'watts' and erdos are options
    N = 64
    '''
    # get adjacency matrix
    conn_mat = load_matrix(conn_method=conn_method)
    print("mat shape:", conn_mat.shape)

    threshold_20_percent_density = 0.137 # 0.04597 for PDC
    adj_mat = compute_adjacency(conn_mat, threshold=threshold_20_percent_density)  
    print("Resulting network density:", np.sum(adj_mat)/4032)

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


    #print('\n================== P 2.3 ==============================')
    #p2_3(freq=10, run='R01', threshold_pdc=0.1226, threshold_dtf=0.1378)

    print('\n================== P 2.4 ==============================')
    p2_4(run="R02")

    #print('\n================== P 2.5 ==============================')
    #G = load_conn_graph(conn="pdc", freq=10, run="R01")
    #p2_5(G)

    #print('\n================== P 2.6 ==============================')
    #p2_6('R01')
    #p2_6('R02')
