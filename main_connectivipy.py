import mne
import connectivipy as cp

PLOTS = False

def save_matrices():
    """
    Save adjacency matrices obtained from DTF and PDC connectivity analysis to file
    """    
    f_dtf = open("data/dtf_matrix.txt", "w")
    f_pdc = open("data/pdc_matrix.txt", "w")

    for i in range(64):
        for j in range(64):
            if j == 63:
                f_dtf.write(str(dtf_significance[i][j]) + "\n")
                f_pdc.write(str(pdc_significance[i][j]) + "\n")
            else:
                f_dtf.write(str(dtf_significance[i][j]) + " ")
                f_pdc.write(str(pdc_significance[i][j]) + " ")

    f_dtf.close()
    f_pdc.close()
    return


file_name = "data/S003R01.edf"
raw_data = mne.io.read_raw_edf(file_name, verbose=True)
print("\nData Info:",raw_data.info)

array_data = raw_data.get_data()
print("array_data shape:", array_data.shape)

# Testing only the first two channels ->['Fc5.', 'Fc3.']
data = cp.Data(array_data, fs=32., chan_names=raw_data.ch_names, data_info='edf_data')
if PLOTS:
    data.plot_data(trial=3)



# fit mvar using Yule-Walker algorithm and order 2,
# you can capture fitted parameters and residual matrix
data.fit_mvar(2, 'yw')
ar, vr = data.mvar_coefficients

#print(cp.conn.conn_estim_dc.keys()) -> connectivity measure available
# investigate connectivity using DTF
dtf_values = data.conn('dtf')
dtf_significance = data.significance(Nrep=200, alpha=0.05)
print("\nPDC sign:",dtf_significance)
if PLOTS:
    data.plot_conn('DTF measure')

# investigate connectivity using PDC
pdc_values = data.conn('pdc')
pdc_significance = data.significance(Nrep=200, alpha=0.05)
print("\nDTF sign:",pdc_significance)
if PLOTS:
    data.plot_conn("PDC measure")

save_matrices()