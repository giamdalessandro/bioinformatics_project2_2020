import warnings
import os
import pyedflib
import numpy as np
import networkx as nx
import connectivipy as cp
import matplotlib.pyplot as plt
import matplotlib.colors as colors


PLOTS        = True
COMPUTE_MATS = False
ADJACENCY    = False
OPTIMIZE_P   = False


def fxn():
    warnings.warn("future",  FutureWarning)
    warnings.warn("warning", Warning)


## point 1.1 func ##

def save_matrices(dtf_mat, pdc_mat, n_channels=64, freq=8, run="R01"):
    """
    Save adjacency matrices obtained from DTF and PDC connectivity analysis to file
        - dtf_mat    : connectivity matrix obtained with DTF measure;
        - pdc_mat    : connectivity matrix obtained with PDC measure;
        - n_channels : number of channels in the data, i.e. the resulting matrices 
                dim (n_channels,n_channels);
        - freq       : frequecy value related to the matrix data;
        - run        : the related run of the experiment, one of {'R01','R02'}.
    """
    dtf_path = "data/dtf_{}_{}hz_auto.txt".format(run,freq)
    pdc_path = "data/pdc_{}_{}hz_auto.txt".format(run,freq)

    print("\nSaving DTF and PDC matrices respectively to {} and {}".format(dtf_path,pdc_path))
    f_dtf = open(dtf_path, "w")
    f_pdc = open(pdc_path, "w")

    for i in range(n_channels):
        for j in range(n_channels):
            if j == 63:
                f_dtf.write(str(dtf_mat[i][j]) + "\n")
                f_pdc.write(str(pdc_mat[i][j]) + "\n")
            else:
                f_dtf.write(str(dtf_mat[i][j]) + " ")
                f_pdc.write(str(pdc_mat[i][j]) + " ")

    f_dtf.close()
    f_pdc.close()
    return


def load_matrix(conn_method="pdc", freq=10, run="R01", auto='auto'):
    """
    Load the adjacency matrix from file
        - conn_method : the method used to compute the connectivity matrix, one of {'dtf','pdc'};
        - freq        : the frequqncy value related to the matrix data;
        - run         : the related run of the experiment, one of {'R01','R02'}.
    """
    mat_file = "data/{}_{}_{}hz_{}.txt".format(conn_method, run, freq, auto) 
    mat_list = []
    print("Loading matrix from '{}' ...".format(mat_file))

    with open(mat_file, "r") as f:
        for row in f.readlines():
            mat_list.append(row.strip().split(" "))
        f.close()

    conn_mat = np.array(mat_list, dtype=np.float32)
    np.fill_diagonal(conn_mat, 0.)
    return conn_mat


def compute_adjacency(conn_mat, threshold=0.1226):    
    """
    Compute binary adjacency matrix from the given connectivity matrix.
        - conn_mat  : the connectivity matrix to be binarified;
        - threshold : each cell of the resulting matrix will be considered 1 if its value
                is greater or equal than 'threshold', 0 otherwise. 
    """
    adj_mat = np.zeros(shape=conn_mat.shape)

    for i in range(conn_mat.shape[0]):
        for j in range(conn_mat.shape[1]):
            if conn_mat[i][j] >= threshold and i != j:
                adj_mat[i][j] = 1

    return adj_mat


def load_conn_graph(conn="pdc", freq=10, run="R01", auto="auto", threshold=0.1226):
    """
    Load the connectivity graph from the related connectivity matrix.
        - conn : the method used to compute the connectivity matrix, one of {'dtf','pdc'};
        - freq : the frequqncy value related to the matrix data;
        - run  : the related run of the experiment, one of {'R01','R02'}.
    """
    print("\nInitializing graph from {}-{}-{}hz matrix ...".format(conn,run,freq))
    adj_mat = compute_adjacency(load_matrix(conn_method=conn,freq=freq,run=run,auto=auto), threshold=threshold)

    G = nx.from_numpy_array(adj_mat,create_using=nx.DiGraph)
    # relabel nodes cause there's no labels in the original EDF file
    # (well, there are, but they are not automatically loaded, so...)
    with open("data/channel_locations.txt") as f:
        mapping = {}
        for line in f:
            l = line.split(sep='        ')
            if l[0] != '\ufeff#':
                mapping.update( { int(l[0]) - 1 : str(l[1])})
    return nx.relabel_nodes(G, mapping)


def already_computed():
    """
    If we have already computed the adj mats for point 1.1 do not compute them again\n
    They are computed for both runs for all freq in the alpha band.
    """
    for run in ['R01', 'R02']:
        for freq in range(8, 14):
            dtf_path = "data/dtf_{}_{}hz_auto.txt".format(run, freq)
            pdc_path = "data/pdc_{}_{}hz_auto.txt".format(run, freq)
            if not os.path.isfile(dtf_path) or not os.path.isfile(pdc_path):
                return False
    return True


def print_adj(conn_method='pdc', freq=10, run='R01', threshold=0.1226, auto='auto'):
    """
    Prints adjacency matrix
    ------------------------
    TODO - - to check this values, prolly wrong
    DTF 10hz R01: threshold of 0.1378 network density -> 0.2006 (20.01%)
    PDC 10hz R01: threshold of 0.1226 network density -> 0.2001 (20.01%)
    DTF 10hz R01: ???
    PDC 10hz R01: ???    
    """

    mat = load_matrix(conn_method=conn_method, freq=freq, run=run, auto=auto)
    plt.matshow(mat)
    plt.title("{} adjacency matrix of run {} @{}Hz".format(conn_method, run, freq))
    plt.colorbar()
    plt.show()

    if threshold is not None:
        mat = compute_adjacency(mat, threshold=threshold)
        density = 100*np.sum(mat)/4032
        plt.matshow(mat)
        plt.title("{} binary adjacency matrix of run {} @{}Hz with density = {:.02f}%".format(conn_method, run, freq, density))
        plt.colorbar()
        plt.show()

 
def p1_1(file_name="data/S003R01_fixed", point='1'):
    
    #### Load EEG data from edf file
    if point == '4':        ### <<<<<<<<<<<<<<<<<<
        if not os.path.isfile(file_name + '_dropped.edf'):
            small_group = ['Fp1.', 'Fp2.', 'F7..', 'F3..', 'Fz..',
                       'F4..', 'F8..', 'T7..', 'C3..', 'Cz..',
                       'C4..', 'T8..', 'P7..', 'P3..', 'Pz..',
                       'P4..', 'P8..', 'O1..', 'O2..']
            pyedflib.highlevel.drop_channels(file_name+".edf", to_keep=small_group)
            # NOTE >>>> first time it computes the edf files it crashes, but the file is there... <<<<<<<<<<<<<<<<<<
        file_name = file_name + '_dropped'

    print("\n[1.{}] >> Analyzing file {}".format(point, file_name))
    f = pyedflib.EdfReader(file_name + ".edf")
    n = f.signals_in_file
    signal_labels = f.getSignalLabels()
    sigbufs = np.zeros((n, f.getNSamples()[0]))
    for i in np.arange(n):
        sigbufs[i, :] = f.readSignal(i)

    print("[1.{}] >> Loaded matrix with shape {}".format(point, sigbufs.shape))
    f.close()

    data = cp.Data(sigbufs, fs=160., chan_names=signal_labels, data_info=file_name)
    with warnings.catch_warnings():             # stupid warning about the plot...
        warnings.simplefilter("ignore")
        fxn()
        data.plot_data(trial=3)

    print("[1.{}] >> Optimizing p...".format(point))
    mv = cp.Mvar
    best_p, crit = mv.order_akaike(sigbufs, p_max=30, method='yw')    
    if PLOTS:
        #### MVAR
        # best_p = 5      # best for both R01 and R02 using 64 channels
        # if point == '4':
        #     if file_name == "data/S003R01_fixed_dropped.edf":
        #         best_p = 13
        #     elif file_name == "data/S003R02_fixed_dropped.edf":
        #         best_p = 14

        plt.plot(1+np.arange(len(crit)), crit, 'g')
        plt.title("Model order estimation")
        plt.xlabel("order(p)")
        plt.ylabel("AIC(p)")
        plt.grid()
        plt.show()

    print("[1.{}] >> Best p = {}".format(point, best_p))
    data.fit_mvar(p=best_p, method='yw')
    
    if point == '4':
        return data

    #### Compute connectivity matrices with DTF and PDC measures
    if not already_computed() or COMPUTE_MATS:
        # investigate connectivity using DTF
        dtf_values = data.conn('dtf',resolution=80)
        dtf_significance = data.significance(Nrep=100, alpha=0.05)
        print("[1.1] >> dtf_shape:",dtf_values.shape)
        print("\n[1.1] >> DTF sign:", dtf_significance)
        if PLOTS:
            data.plot_conn("DTF measure")

        # investigate connectivity using PDC
        pdc_values = data.conn('pdc',resolution=80)
        pdc_significance = data.significance(Nrep=100, alpha=0.05)
        print("[1.1] >> pdc_shape:", pdc_values.shape)
        print("\n[1.1] >> PDC sign:", pdc_significance)
        if PLOTS:
            data.plot_conn("PDC measure")

        for i in range(8,14):
            save_matrices(dtf_mat=dtf_values[i],pdc_mat=pdc_values[i],n_channels=64,freq=i,run=file_name[9:12])

    if PLOTS:
        print_adj(conn_method='pdc', freq=10, run='R01', threshold=0.1226)
        print_adj(conn_method='dtf', freq=10, run='R01', threshold=0.1378)


def p1_2():
    print("[1.2] >> Still to be implemented...")


def p1_3():
    print("[1.3] >> Still to be implemented...")


def p1_4(R='R01'):
    data = p1_1(point='4')
    pdc_path = "data/pdc_{}_10hz_19_channels.txt".format(R)

    if not os.path.isfile(pdc_path):
        pdc_values = data.conn('pdc', resolution=80)
        data.significance(Nrep=100, alpha=0.05)         # returns pdc_significance but what to do with it? However, it does side effect
        pdc_mat = pdc_values[10]
        f_pdc = open(pdc_path, "w")
        for i in range(19):
            for j in range(19):
                if j == 18:
                    f_pdc.write(str(pdc_mat[i][j]) + "\n")
                else:
                    f_pdc.write(str(pdc_mat[i][j]) + " ")
        f_pdc.close()
    
    if PLOTS:
        print("[1.4] >> Plotting connectivity...")
        print_adj(conn_method='pdc', freq=10, run='R01', threshold=None, auto='19_channels')
        # data.plot_conn("PDC measure")     # is it even useful?


def load_channel_coordinates(label=True, map_ch=False):
    """
    Loads channels coordinates in a disctionary and returns it
    """
    with open("data/channel_locations.txt") as f:
        pos = {}
        for line in f:
            # yes, there are 8 spaces in the file.
            l = line.split(sep='        ')
            if l[0] != '\ufeff#':
                if label:
                    pos.update({str(l[1]): [float(l[2]), float(l[3])]})
                elif map_ch:
                    pos.update({int(l[0])-1 : str(l[1])})
                else:
                    pos.update({int(l[0])-1 : [float(l[2]), float(l[3])]})
    return pos


def p1_5(G, point='1.5', communities=None, nodelist=None, edgelist=None):
    """
    Prints a topological representation of the networks
    Node colors depend on their degree
    """
    if nodelist is None:
        nodelist = G.nodes()
    if edgelist is None:
        edgelist = G.edges()
    
    pos = load_channel_coordinates()

    def p1_5_helper(G, pos, degree):
        """
        Helper function to now write two times the same plt stuff
        """
        node_color = []
        for node in G.nodes():
            if degree == 'in':
                node_color.append(G.in_degree(node))
            else:
                node_color.append(G.out_degree(node))

        cmap = 'viridis'
        vmin = min(node_color)
        vmax = max(node_color)

        nx.draw_networkx(G, pos=pos, arrows=True, with_labels=True, vmin=vmin, vmax=vmax,
                         node_size=700, edge_color='black', node_color=node_color, cmap=cmap,
                         edgelist=edgelist, nodelist=nodelist)

        plt.title("Topological representation of the network - {} degree".format(degree))

        sm = plt.cm.ScalarMappable(
            cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
        sm._A = []
        plt.colorbar(sm)
        plt.show()


    def p4_3_helper(G, pos, communities):

        cmap = colors.ListedColormap(
            ['#1f78b4', '#33a02c', '#e31a1c', '#ff7f00', '#6a3d9a'], 'indexed', max(communities) + 1)
        # vmin = min(communities)
        # vmax = max(communities)
        nx.draw_networkx(G, pos=pos, arrows=True, with_labels=True,# vmin=vmin, vmax=vmax,
                         node_size=700, edge_color='black', node_color=communities, cmap=cmap)
        plt.title("Topological representation of the network's communities found with Infomap ")
        plt.show()


    if point == '1.5':
        p1_5_helper(G, pos, 'in')
        p1_5_helper(G, pos, 'out')
    elif point == '4.3':
        p4_3_helper(G, pos, communities)

def p1_6():
    print("[1.6] >> Still to be implemented...")


### MAIN 
#p1_4()
