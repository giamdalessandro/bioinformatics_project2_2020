import os
import pyedflib
import numpy as np
import connectivipy as cp
import matplotlib.pyplot as plt

'''
p = []
#### Loading EEG data
for f in os.listdir("test/"):
    file_name = "test/" + f
    print("\nAnalyzing file", file_name)
    raw_data = mne.io.read_raw_edf(file_name, verbose=True)
    print("Data Info:",raw_data.info)

    events, event_dict = mne.events_from_annotations(raw_data)
    #print("Events dict:",event_dict)
    #print("Read events:",events)

    array_data = raw_data.get_data()
    print("array_data shape:", array_data.shape)

    #### Model order
    mv = cp.Mvar
    # find best model order using Vieira-Morf algorithm
    best_p, crit = mv.order_akaike(array_data, p_max=15, method='yw')
    #best_p = 12
    #plt.plot(1+np.arange(len(crit)), crit, 'g')
    #plt.title("Model order estimation")
    #plt.xlabel("order(p)")
    #plt.ylabel("AIC(p)")
    #plt.grid()
    #plt.show()
    print(crit)
    print(best_p)

    p.append(best_p)

print(p)
'''



file_name = "test/S050R02.edf"
print("\nAnalyzing file", file_name)
print(pyedflib.highlevel.read_edf_header(file_name, read_annotations=True))

f = pyedflib.EdfReader(file_name)
#n = f.open()
#print("number of signals in file:", n)
#signal_labels = f.getSignalLabels()
#print("signal_labels:", signal_labels)

f.close()
