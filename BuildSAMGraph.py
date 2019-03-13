# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 11:33:00 2019

Build SAM script is a precursor to Monte Carlo Simulation
This script is developed for a full/continuous SAM surface

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

## GRID SETUP: --------------------------------------------------------

# Lattice Parameters:
lattice_a = 0.5
lattice_b = 0.25
lattice_c = 0.433

x_size = 10.4
y_size = 11*0.9

fig, ax = plt.subplots(figsize=(5,4))

# Setup Grid:
dx,dy = lattice_a,lattice_c
x,y = np.mgrid[slice(0,x_size+dx,dx),slice(0,y_size+dy,dy)]
if np.remainder(np.shape(x)[0],2)!=0:
    print('INCREASE x BY dx',np.remainder(np.shape(x)[0],2))
if np.remainder(np.shape(y)[0],2)!=0:
    print('INCREASE y BY dy',np.remainder(np.shape(y)[0],2))


# Add Offset:
for i in range(0,np.shape(x)[1]):
    index = np.remainder(i,2)
    if index == 1:
        x[:,i] = x[:,i]+ lattice_b

# Add Ligand ID:
strand_id = x*0
for i in range(0,np.shape(strand_id)[0]):
    for j in range(0,np.shape(strand_id)[1]):
        index1 = np.remainder(j,2) #+np.remainder(j,4)
        index2 = 0
        if np.remainder(j,4) == 1:
            index2 = np.remainder(i,2)
        if np.remainder(j,4) == 3:
            index2 = np.remainder(i,2)-1            
        if (index1 == 1)&(index2 == 0):
            strand_id[i,j] = 1
# Plot:
ax.scatter(x,y,c=strand_id,s=10)
k=0
for i in range(0,np.shape(strand_id)[0]):
    for j in range(0,np.shape(strand_id)[1]):
        if strand_id[i,j]==1:
            k = k+1
            ax.annotate(k-1,(x[i,j],y[i,j]))


plt.show()

## GRAPH SETUP: ------------------------------------------------------------

# Calculate the number of ligands in the system:
nlig = int(np.sum(strand_id))
print 'Number of Ligands is', nlig

class Ligand:
    def __init__(self, x_coord,y_coord,lig_index):
        self.x = x_coord
        self.y = y_coord
        self.index = lig_index

lig_index = 0
ligands = []
for i in range(0,np.shape(strand_id)[0]):
    for j in range(0,np.shape(strand_id)[1]):
        if strand_id[i,j] == 1:
            single_lig = Ligand(x[i,j],y[i,j],lig_index) 
            lig_index = lig_index + 1
            ligands.append(single_lig)
            
x_size_true = x_size +0.6
y_size_true = y_size+0.5

sam_full = Graph()
sam_full.add_vertices(int(np.sum(strand_id)))

for i in range(nlig-1):
    for j in range(i+1,nlig):
        xdist = np.abs(ligands[i].x-ligands[j].x)
        ydist = np.abs(ligands[i].y-ligands[j].y)

        # Account for periodic boundary conditions:
        if xdist > 0.5*x_size_true:
            xdist = x_size_true-xdist
        if ydist > 0.5*y_size_true:
            ydist = y_size_true-ydist
        dist = np.sqrt(xdist**2+ydist**2)
        
        if abs(dist-1)<0.01:
            sam_full.add_edges([(i,j)])


