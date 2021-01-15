import os
import warnings
import pyedflib
import numpy as np
import networkx as nx
import connectivipy as cp
from itertools import count
import matplotlib.pyplot as plt
import matplotlib.colors as colors

PLOTS        = True
COMPUTE_MATS = False
ADJACENCY    = False
OPTIMIZE_P   = False



THRES_DTF_10HZ_R01_20percent = 0.13785000000003989
THRES_DTF_10HZ_R02_20percent = 0.1323000000000405

THRES_PDC_25HZ_R01_20percent = 0.1268000000000301
THRES_PDC_25HZ_R02_20percent = 0.12200000000003064

THRES_PDC_10HZ_R01_01percent = 0.26535000000000053
THRES_PDC_10HZ_R01_05percent = 0.18070000000000985
THRES_PDC_10HZ_R01_10percent = 0.152150000000013
THRES_PDC_10HZ_R01_20percent = 0.12270000000001624
THRES_PDC_10HZ_R01_30percent = 0.09830000000001893
THRES_PDC_10HZ_R01_50percent = 0.06645000000002244

THRES_PDC_10HZ_R02_01percent = 0.261250000000001
THRES_PDC_10HZ_R02_05percent = 0.17730000000001023
THRES_PDC_10HZ_R02_10percent = 0.14730000000001353
THRES_PDC_10HZ_R02_20percent = 0.1167500000000169
THRES_PDC_10HZ_R02_30percent = 0.09535000000001925
THRES_PDC_10HZ_R02_50percent = 0.061450000000022986

RUNS = ['R01', 'R02']


def fxn():
    warnings.warn("future",  FutureWarning)
    warnings.warn("warning", Warning)
