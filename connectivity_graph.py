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
        - dtf_mat   : connectivity matrix obtained with DTF measure;
        - pdc_mat   : connectivity matrix obtained with PDC measure;
        - n_channels: number of channels in the data, i.e. the resulting matrices 
                dim (n_channels,n_channels);
        - freq      : frequecy value of analysis.
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

def load_matrix(conn_method="dtf", freq=10, run="R01"):
    """
    Load the adjacency matrix from file
        - conn_method: the method used to compute the connectivity matrix, one of {'dtf','pdc'};
        - freq       : the frequqncy value related to the matrix data;
        - run        : the related run of the experiment, one of {'R01','R02'}.
    """
    mat_file = "data/{}_{}_{}hz_auto.txt".format(conn_method,run,freq) 
    mat_list = []
    print("Loading matrix from '{}' ...".format(mat_file))

    with open(mat_file, "r") as f:
        for row in f.readlines():
            mat_list.append(row.strip().split(" "))
        f.close()

    return np.array(mat_list, dtype=np.float32)

def compute_adjacency(conn_mat, threshold=0.05):    
    """
    Compute binary adjacency matrix from the given connectivity matrix.
        - conn_mat : the connectivity matrix to be binarified;
        - threshold: each cell of the resulting matrix will be considered 1 if its value
                is greater or equal than 'threshold', 0 otherwise. 
    """
    adj_mat = np.zeros(shape=conn_mat.shape)

    for i in range(conn_mat.shape[0]):
        for j in range(conn_mat.shape[1]):
            if conn_mat[i][j] >= threshold and i != j:
                adj_mat[i][j] = 1

    return adj_mat

def load_conn_graph(conn="dtf", freq=10, run="R01"):
    """
    Load the connectivity graph from the related connectivity matrix.
        - conn : the method used to compute the connectivity matrix, one of {'dtf','pdc'};
        - freq : the frequqncy value related to the matrix data;
        - run  : the related run of the experiment, one of {'R01','R02'}.
    """
    print("\nInitializing graph from {}-{}-{}hz matrix ...".format(conn,run,freq))
    adj_mat = compute_adjacency(load_matrix(conn_method=conn,freq=freq,run=run))

    return nx.from_numpy_array(adj_mat,create_using=nx.DiGraph)



file_name = "data/S003R02_fixed.edf"
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
# find best model order using Vieira-Morf algorithm
best_p, crit = mv.order_akaike(sigbufs, p_max=15, method='yw')
if PLOTS:
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

    for i in range(8,14):
        save_matrices(dtf_mat=dtf_values[i],pdc_mat=pdc_values[i],n_channels=64,freq=i,run=file_name[9:12])


#### Compute adjacency matrix 
"""
TODO -- to check this values, prolly wrong
DTF: with a threshold of 0.07881 we obtain a neetwork density of 0.2006 (20.01%) 
PDC: with a threshold of 0.04597 we obtain a neetwork density of 0.2003 (20.03%)
"""
if ADJACENCY:
    conn_mat = load_matrix(conn_method="dtf")
    print("mat shape:",conn_mat.shape)

    adj_mat = compute_adjacency(conn_mat, threshold=0.03)  # 0.04597 for PDC
    #print(adj_mat)
    max_edges = 64*(64-1)
    print("Resutling network density:", np.sum(adj_mat)/max_edges)


"""
print()
print(dtf_values[10])
print()
print(pdc_values[10])
"""
