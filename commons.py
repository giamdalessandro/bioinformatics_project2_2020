import warnings
import os
import pyedflib
import numpy as np
import networkx as nx
import connectivipy as cp
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from itertools import count

PLOTS        = True
COMPUTE_MATS = False
ADJACENCY    = False
OPTIMIZE_P   = False


THRES_10HZ_R01_20percent = 0.12270000000001624
THRES_10HZ_R02_20percent = 0.1167500000000169


def fxn():
    warnings.warn("future",  FutureWarning)
    warnings.warn("warning", Warning)
