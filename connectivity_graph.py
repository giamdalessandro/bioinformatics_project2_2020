import os
import pyedflib
import numpy as np
import networkx as nx
import connectivipy as cp
import matplotlib.pyplot as plt

PLOTS        = False
COMPUTE_MATS = False
ADJACENCY    = False

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

def load_matrix(conn_method="pdc", freq=10, run="R01"):
    """
    Load the adjacency matrix from file
        - conn_method : the method used to compute the connectivity matrix, one of {'dtf','pdc'};
        - freq        : the frequqncy value related to the matrix data;
        - run         : the related run of the experiment, one of {'R01','R02'}.
    """
    mat_file = "data/{}_{}_{}hz_auto.txt".format(conn_method,run,freq) 
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

def load_conn_graph(conn="pdc", freq=10, run="R01"):
    """
    Load the connectivity graph from the related connectivity matrix.
        - conn : the method used to compute the connectivity matrix, one of {'dtf','pdc'};
        - freq : the frequqncy value related to the matrix data;
        - run  : the related run of the experiment, one of {'R01','R02'}.
    """
    print("\nInitializing graph from {}-{}-{}hz matrix ...".format(conn,run,freq))
    adj_mat = compute_adjacency(load_matrix(conn_method=conn,freq=freq,run=run))

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


def p1_1_print_adj():
    """
    Prints adjacency matrix
    """

    ## PDC ##
    mat = load_matrix(conn_method='pdc', freq=10, run='R01')
    plt.matshow(mat)
    plt.title("PDC adjacency matrix")
    plt.colorbar()
    plt.show()

    mat = compute_adjacency(mat, threshold=0.1226)
    density = 100*np.sum(mat)/4032
    plt.matshow(mat)
    plt.title("PDC binary adjacency matrix with density = {:.02f}%".format(density))
    plt.colorbar()
    plt.show()

    ## DTF ##
    mat = load_matrix(conn_method='dtf', freq=10, run='R01')
    plt.matshow(mat)
    plt.title("DTF adjacency matrix")
    plt.colorbar() 
    plt.show()

    mat = compute_adjacency(mat, threshold=0.1378)
    density = 100*np.sum(mat)/4032
    plt.matshow(mat)
    plt.title("DTF binary adjacency matrix with density = {:.02f}%".format(density))
    plt.colorbar()
    plt.show()

 
def p1_1(file_name="data/S003R02_fixed", point='1'):
    
    #### Load EEG data from edf file
    if point == '4':        ### <<<<<<<<<<<<<<<<<<
        if not os.path.isfile(file_name + '_dropped.edf'):
            small_group = ['Fp1.', 'Fp2.', 'F7..', 'F3..', 'Fz..',
                       'F4..', 'F8..', 'T7..', 'C3..', 'Cz..',
                       'C4..', 'T8..', 'P7..', 'P3..', 'Pz..',
                       'P4..', 'P8..', 'O1..', 'O2..']
            pyedflib.highlevel.drop_channels(file_name+".edf", to_keep=small_group)
            # first time it computes the edf files it crashes, but the file is there... <<<<<<<<<<<<<<<<<<
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
    data.plot_data(trial=3)

    #### MVAR
    if point == '1':    # best for both R01 and R02
        best_p = 5
    elif point == '4':
        if file_name == "data/S003R01_fixed_dropped.edf":
            best_p = 13
        elif file_name == "data/S003R02_fixed_dropped.edf":
            best_p = 14
    
    if PLOTS:
        mv = cp.Mvar
        best_p, crit = mv.order_akaike(sigbufs, p_max=30, method='yw')
        plt.plot(1+np.arange(len(crit)), crit, 'g')
        plt.title("Model order estimation")
        plt.xlabel("order(p)")
        plt.ylabel("AIC(p)")
        plt.grid()
        plt.show()

    print("[1.{}] >> Best p = {}".format(point, best_p))
    data.fit_mvar(p=best_p, method='yw')
    #if PLOTS:
    p1_1_print_adj()

    if point == '4':
        return data

    #### Compute connectivity matrices with DTF and PDC measures
    if COMPUTE_MATS:
        # investigate connectivity using DTF
        dtf_values = data.conn('dtf',resolution=80)
        dtf_significance = data.significance(Nrep=100, alpha=0.05)
        print("dtf_shape:",dtf_values.shape)
        print("\nDTF sign:",dtf_significance)
        if PLOTS:
            data.plot_conn("DTF measure")

        # investigate connectivity using PDC
        pdc_values = data.conn('pdc',resolution=80)
        #pdc_values = data.short_time_conn('pdc', nfft=100, no=10)
        pdc_significance = data.significance(Nrep=100, alpha=0.05)
        print("pdc_shape:",pdc_values.shape)
        print("\nPDC sign:",pdc_significance)
        if PLOTS:
            data.plot_conn("PDC measure")
            #data.plot_short_time_conn("PDC")

        #for i in range(8,14):
        #    save_matrices(dtf_mat=dtf_values[i],pdc_mat=pdc_values[i],n_channels=64,freq=i,run=file_name[9:12])

    """
    TODO -- to check this values, prolly wrong
    DTF 10hz R01: threshold of 0.1378 network density -> 0.2006 (20.01%) 
    PDC 10hz R01: threshold of 0.1226 network density -> 0.2001 (20.01%)
    """


def p1_2():
    print("[1.2] >> Still to be implemented...")


def p1_3():
    print("[1.3] >> Still to be implemented...")


def p1_4(R='R01'):
    data = p1_1(point='4')
    pdc_values = data.conn('pdc', resolution=80)
    pdc_significance = data.significance(Nrep=100, alpha=0.05)
    pdc_mat = pdc_values[10]
    pdc_path = "data/pdc_{}_10hz_19_channels.txt".format(R)

    if not os.path.isfile(pdc_path):
        f_pdc = open(pdc_path, "w")
        for i in range(19):
            for j in range(19):
                if j == 18:
                    f_pdc.write(str(pdc_mat[i][j]) + "\n")
                else:
                    f_pdc.write(str(pdc_mat[i][j]) + " ")
        f_pdc.close()
    
    if PLOTS:
        data.plot_conn("PDC measure")


def p1_5(G, nodelist=None, edgelist=None):
    """
    Prints a topological representation of the networks
    Node colors depend on their degree
    """
    if nodelist is None:
        nodelist = G.nodes()
    if edgelist is None:
        edgelist = G.edges()
    with open("data/channel_locations.txt") as f:
        pos = {}
        for line in f:
            # yes, there are 8 spaces in the file.
            l = line.split(sep='        ')
            if l[0] != '\ufeff#':
                pos.update({str(l[1]): [float(l[2]), float(l[3])]})

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

        plt.title(
            "Topological representation of the network - {} degree".format(degree))

        sm = plt.cm.ScalarMappable(
            cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
        sm._A = []
        plt.colorbar(sm)
        plt.show()

    p1_5_helper(G, pos, 'in')
    p1_5_helper(G, pos, 'out')



### MAIN 
p1_1()
#p1_4()
