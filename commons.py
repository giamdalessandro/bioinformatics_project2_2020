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


THRES_10HZ_R01_20percent = 0.12270000000001624
THRES_10HZ_R02_20percent = 0.1167500000000169


RUNS = ['R01', 'R02']


def fxn():
    warnings.warn("future",  FutureWarning)
    warnings.warn("warning", Warning)
