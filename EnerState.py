# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:46:56 2019

@author: camil
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import matplotlib.colors as col
from scipy.optimize import curve_fit
import pandas as pd
from math import factorial
import igraph
from igraph import *
import random
import time


## DEFINE ENERGY FUNCTION: --------------------------------------------------

def EnerState(graph,A,B):
    w1 = float(np.shape(graph.get_edgelist())[0])
    w2 = np.shape(graph.maximal_cliques(min=3,max=3))[0]
    ener = w1*A-w2*B  #w1*0.25-w2*1.5
    return ener