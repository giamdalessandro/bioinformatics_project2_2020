import mne
from pdc_tdf import *

file_name = "S003R01.edf"
data = mne.io.read_raw_edf(file_name, verbose=True)

# you can get the metadata included in the file_name and a list of all channels:
raw_data = data.get_data()
info = data.info
channels = data.ch_names
print("\n",info)
print("\n",channels)

#data.plot_psd(fmax=50)
data.plot(start=0,duration=5, n_channels=30, verbose=True)

"""
#events = mne.read_events(file_name + ".event")
#print(events)
"""