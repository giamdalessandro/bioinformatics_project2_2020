import mne
import numpy as np
import connectivipy as cp

PLOTS        = True
COMPUTE_MATS = False

def save_matrices(dtf_mat, pdc_mat, n_channels=64):
    """
    Save adjacency matrices obtained from DTF and PDC connectivity analysis to file
        - n_channels: number of channels in the data, i.e. the resulting matrices 
                dim (n_channels,n_channels).
    """    
    print("\nSaving DTF and PDC matrices respectively to 'data/dtf_matrix.txt' and 'data/pdc_matrix.txt'")
    f_dtf = open("data/dtf_matrix.txt", "w")
    f_pdc = open("data/pdc_matrix.txt", "w")

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

def load_matrix(conn_method="DTF"):
    """
    Load the adjacency matrix from file
        - conn_method: the method used to compute the connectivity matrix, one of {'DTF','PDC'}.
    """
    mat_file = "data/dtf_matrix.txt" if conn_method == "DTF" else "data/pdc_matrix.txt" 
    mat_list = []
    print("\nLoading matrix from '{}' ...".format(mat_file))

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
                is >= than threshold, 0 otherwise. 
    """
    adj_mat = np.zeros(shape=conn_mat.shape)

    for i in range(conn_mat.shape[0]):
        for j in range(conn_mat.shape[1]):
            if conn_mat[i][j] >= threshold and i != j:
                adj_mat[i][j] = 1

    return adj_mat

#### Loading EEG data
file_name = "data/S003R01.edf"
raw_data = mne.io.read_raw_edf(file_name, verbose=True)
print("\nData Info:",raw_data.info)

array_data = raw_data.get_data()
print("array_data shape:", array_data.shape)

data = cp.Data(array_data, fs=160., chan_names=raw_data.ch_names, data_info='edf_data')
if PLOTS:
    data.plot_data(trial=3)

#### Compute connectivity matrices with DTF and PDC
if COMPUTE_MATS:
    # fit mvar using Yule-Walker algorithm and order 2,
    # you can capture fitted parameters and residual matrix
    data.fit_mvar(2, 'yw')
    ar, vr = data.mvar_coefficients

    # investigate connectivity using DTF
    dtf_values = data.conn('dtf')
    dtf_significance = data.significance(Nrep=100, alpha=0.05)
    print(dtf_values[0])
    print("\nDTF sign:",dtf_significance)
    if PLOTS:
        data.plot_conn('DTF measure')

    # investigate connectivity using PDC
    pdc_values = data.conn('pdc')
    pdc_significance = data.significance(Nrep=100, alpha=0.05)
    print(pdc_values[0])
    print("\nPDC sign:",pdc_significance)
    if PLOTS:
        data.plot_conn("PDC measure")

    save_matrices(dtf_significance,pdc_significance,n_channels=64)


#### Compute adjacency matrix 
"""
DTF: with a threshold of 0.07881 we obtain a neetwork density of 0.2006 (20.01%) 
PDC: with a threshold of 0.04597 we obtain a neetwork density of 0.2003 (20.03%)

conn_mat = load_matrix(conn_method='DTF')
print("mat shape:",conn_mat.shape)

adj_mat = compute_adjacency(conn_mat, threshold=0.07881)  # 0.04597 for PDC
#print(adj_mat)
max_edges = 64*(64-1)
print("Network density:", np.sum(adj_mat)/max_edges)
"""