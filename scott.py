import mne
import scot
import numpy as np
import connectivipy as cp
import matplotlib.pyplot as plt

PLOTS = False
COMPUTE_MATS = True
ADJACENCY = True

#### Loading EEG data
file_name = "data/S003R01.edf"
print("\nAnalyzing file", file_name)
raw_data = mne.io.read_raw_edf(file_name, verbose=True)
print("Data Info:", raw_data.info)

events, event_dict = mne.events_from_annotations(raw_data)
print("Events dict:", event_dict)
print("Read events:", events)

array_data = raw_data.get_data()
print("array_data shape:", array_data.shape)


mvar = scot.var.VAR(1)
mvar.fit(array_data)
print(mvar.p)
A = [array_data[:, :2440],
     array_data[:, 2440:4880],
     array_data[:, 4880:7320],
     array_data[:, 7320:]]

mvar.optimize_order(
    data=A, min_p=1, max_p=15, n_jobs=-1, verbose=2)
print(mvar.p)
# data : array-like, shape (n_trials, n_samples, n_channels)
#        Segmented data set on which to optimize the model order. At least 2 trials are required
