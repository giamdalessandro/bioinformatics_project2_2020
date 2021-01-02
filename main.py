import mne
from scot.var import VAR
from scot.connectivity import Connectivity

file_name = "data/S003R01.edf"
data = mne.io.read_raw_edf(file_name, verbose=True)

# you can get the metadata included in the file_name and a list of all channels:
info = data.info
channels = data.ch_names
print("\n",info)
print("\n",channels)

raw_data = data.get_data()
print("raw data shape :", raw_data.shape)

"""
#data.plot_psd(fmax=50)
#data.plot(start=0,duration=5, n_channels=30, verbose=True)
"""

var = VAR(model_order=40)
fitted = var.fit(raw_data)
print("VAR coefficient:",fitted.coef.shape)

conn = Connectivity(fitted.coef)