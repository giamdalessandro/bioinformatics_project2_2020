import mne
import matplotlib.pyplot as plt

from scot.var import VAR
from scot.connectivity import Connectivity
from scot.plotting import plot_connectivity_spectrum, plot_connectivity_significance
from scot.eegtopo.topoplot import Topoplot

"""
To check also connectiviPy module, maybe better
"""

file_name = "data/S003R01.edf"
data = mne.io.read_raw_edf(file_name, verbose=True)

# get the metadata included in the file_name and a list of all channels:
info = data.info
channels = data.ch_names
print("\n",info)
print("\n",channels)

raw_data = data.get_data()
print("raw_data shape:", raw_data.shape)

"""
#data.plot_psd(fmax=50)
#data.plot(start=0,duration=5, n_channels=30, verbose=True)
"""

var = VAR(model_order=40)
fitted = var.fit(raw_data)
print("\nVAR coefficient:",fitted.coef.shape)

conn = Connectivity(fitted.coef)
pdc  = conn.PDC()
print("PDC shape:", pdc.shape)


#print("Channel {} topoplot...".format(channels[0]))
#plot_circular(pdc[:,:,0], colors=(0,255,0), topo=Topoplot(), topomaps=[80,80])

# plotting only 3 channel
print("Plotting {} channels...".format(channels[:3]))
plot_connectivity_spectrum(pdc[:3,:3,:], freq_range=(0,80))
#plot_connectivity_significance(pdc[:3,:3,:], freq_range=(0,80))
plt.show()