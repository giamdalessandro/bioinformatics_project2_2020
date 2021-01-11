import pyedflib
import numpy as np
import networkx as nx
import connectivipy as cp
import matplotlib.pyplot as plt

PLOTS        = False
COMPUTE_MATS = False
ADJACENCY    = False


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

    return np.array(mat_list, dtype=np.float32)

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


def print_adj():
    """
    Prints adjacency matrix
    """

    ## PDC ##
    mat = load_matrix(conn_method='pdc', freq=10, run='R01')
    plt.matshow(mat)
    plt.title("PDC adjacency matrix")
    plt.show()

    mat = compute_adjacency(mat, threshold=0.1226)
    density = 100*np.sum(mat)/4032
    plt.matshow(mat)
    plt.title("PDC binary adjacency matrix with density = {:.02f}%".format(density))
    plt.show()

    ## DTF ##
    mat = load_matrix(conn_method='dtf', freq=10, run='R01')
    plt.matshow(mat)
    plt.title("DTF adjacency matrix")
    plt.show()

    mat = compute_adjacency(mat, threshold=0.1378)
    density = 100*np.sum(mat)/4032
    plt.matshow(mat)
    plt.title("DTF binary adjacency matrix with density = {:.02f}%".format(density))
    plt.show()



def p1_5(G):
    with open("data/channel_locations.txt") as f:
        pos = {}
        for line in f:
            l = line.split(sep='        ')
            if l[0] != '\ufeff#':
                pos.update({ str(l[1]) : [float(l[2]), float(l[3])] })

    nx.draw_networkx(G, pos=pos, arrows=True, with_labels=True)




if __name__ == "__main__":
    #### Loading EEG data from edf file
    file_name = "data/S003R01_fixed.edf"
    print("\nAnalyzing file", file_name)
    f = pyedflib.EdfReader(file_name)
    n = f.signals_in_file
    signal_labels = f.getSignalLabels()
    sigbufs = np.zeros((n, f.getNSamples()[0]))
    for i in np.arange(n):
        sigbufs[i, :] = f.readSignal(i)

    print("Loaded matrix with shape", sigbufs.shape)
    f.close()


    data = cp.Data(sigbufs, fs=160., chan_names=signal_labels, data_info=file_name)
    if PLOTS:
        data.plot_data(trial=3)


    #### Model order
    mv = cp.Mvar
    best_p = 5
    if PLOTS:
        best_p, crit = mv.order_akaike(sigbufs, p_max=15, method='yw')
        plt.plot(1+np.arange(len(crit)), crit, 'g')
        plt.title("Model order estimation")
        plt.xlabel("order(p)")
        plt.ylabel("AIC(p)")
        plt.grid()
        plt.show()
        print(crit)

    print("Best p =", best_p)


    # fit mvar using Yule-Walker algorithm and order p
    data.fit_mvar(p=best_p, method='yw')
    ar, vr = data.mvar_coefficients
    print(data._parameters)


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

    G = load_conn_graph(conn="pdc", freq=10, run="R01")
    p1_5(G)
