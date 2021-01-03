import mne
import connectivipy as cp


file_name = "data/S003R01.edf"
raw_data = mne.io.read_raw_edf(file_name, verbose=True)
print("\nData Info:",raw_data.info)

array_data = raw_data.get_data()
print("array_data shape:", array_data.shape)

# Testing only the first two channels
data = cp.Data(array_data[:2,:], fs=32., chan_names=['Fc5.', 'Fc3.'], data_info='edf_data')
data.plot_data(trial=3)



# fit mvar using Yule-Walker algorithm and order 2
data.fit_mvar(2, 'yw')

# you can capture fitted parameters and residual matrix
ar, vr = data.mvar_coefficients

#print(cp.conn.conn_estim_dc.keys()) -> connectivity measure available
# investigate connectivity using DTF
dtf_values = data.conn('dtf')
dtf_significance = data.significance(Nrep=200, alpha=0.05)
data.plot_conn('DTF measure')
print("\nPDC sign:",dtf_significance)

# investigate connectivity using PDC
pdc_values = data.conn('pdc')
pdc_significance = data.significance(Nrep=200, alpha=0.05)
data.plot_conn("PDC measure")
print("\nDTF sign:",pdc_significance)