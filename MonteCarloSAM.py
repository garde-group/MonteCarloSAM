# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 16:48:54 2019

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

def MonteCarloSAM(graph, nsteps, A, B):
    
    ## RUN MONTE CARLO SIMULATION: -----------------------------------------------
    
    start = time.time()
    
    # Monte Carlo Simulation:
    
    kT = 2.48 # Gives units of kJ/mol
    
    # Create simulated graph
    edge_list = graph.get_edgelist()
    sim_graph_0 = Graph()
    full_EL = graph.get_edgelist() # Full list of possible edges to make or break
    
    # Initialize monitoring variables:
    prob_mon = np.zeros(nsteps)
    ener_0_mon = np.zeros(nsteps)
    ener_diff_mon = np.zeros(nsteps)
    move_mon = np.zeros(nsteps)
    break_mon = np.zeros(nsteps)  
    n_edges_mon = np.zeros(nsteps)
    clique_mon = np.zeros((3,nsteps))
    p_null = np.zeros(nsteps)
    p_dimer = np.zeros(nsteps)
    p_double_dimer = np.zeros(nsteps)
    p_trimer1 = np.zeros(nsteps)
    p_trimer2 = np.zeros(nsteps)
    p_tetramer = np.zeros(nsteps)
    
    for i in range(0,nsteps):  
        if np.mod(i,1000)==0:
            print('Monte Carlo Step', i)
            
        # Create simulated edge list:
        sim_EL_0 = sim_graph_0.get_edgelist()
        
        # Select random edge:
        edge_sel = random.randint(0,np.shape(full_EL)[0]-1)
        
        # Set initial energy:
        ener_0 = EnerState(sim_graph_0,A,B)
        
        # Calculate graph after move:
        # Does selected edge exist in simulated graph?:
    
        if not any([full_EL[edge_sel]==sim_EL_0[j] for j in range(0,np.shape(sim_EL_0)[0])]):
            # If not, add this edge:
            sim_EL_1 = list(sim_EL_0)
            sim_EL_1.append(full_EL[edge_sel])
        elif any([full_EL[edge_sel]==sim_EL_0[j] for j in range(0,np.shape(sim_EL_0)[0])]):
            # If so, remove this edge:
            sim_EL_1 = list(sim_EL_0)
            sim_EL_1.remove(full_EL[edge_sel])
        
        sim_graph_1 = Graph(sim_EL_1)
            
        ener_1 = EnerState(sim_graph_1,A,B)
        
        # Calculate the probability of making the move:
        prob=np.exp(-(ener_1-ener_0)/kT)/(1+np.exp(-(ener_1-ener_0)/kT))
        
        criterion = random.uniform(0, 1)
        if criterion < prob:
            sim_graph_0 = Graph(sim_EL_1)
            
        
        # Monitoring variables:
        prob_mon[i] = prob
        ener_0_mon[i] = ener_0
        ener_diff_mon[i] = ener_1-ener_0
        move_mon[i] = criterion < prob
        break_mon[i] = any([full_EL[edge_sel]==sim_EL_0[j] for j in range(0,np.shape(sim_EL_0)[0])])
        n_edges_mon[i] = np.shape(sim_EL_0)[0]
        clique_mon[0,i] = np.shape(sim_graph_0.maximal_cliques(min=1,max=1))[0]
        clique_mon[1,i] = np.shape(sim_graph_0.maximal_cliques(min=2,max=2))[0]
        clique_mon[2,i] = np.shape(sim_graph_0.maximal_cliques(min=3,max=3))[0]
        np.shape(edge_list)[0]
        
        # If graph is small (non-continuous):
        if np.shape(edge_list)[0]<10:
            p_null[i] = int(n_edges_mon[i]==0)
            p_dimer[i] = int(n_edges_mon[i]==1)
            p_double_dimer[i] = int(n_edges_mon[i]==2)
            p_trimer1[i] = 0
            p_trimer2[i] = 0
            if np.any([j==(0,3) for j in sim_EL_1]):
                if np.any([j==(0,2) for j in sim_EL_1])&np.any([j==(2,3) for j in sim_EL_1]):
                    p_trimer1[i] = 1
                if np.any([j==(0,1) for j in sim_EL_1])&np.any([j==(1,3) for j in sim_EL_1]):
                    p_trimer2[i] = 1
                    
                p_tetramer[i] = int(n_edges_mon[i]==5)
            
        
    end = time.time()
    
    print 'Time Elapsed:', end-start
    
    return clique_mon,n_edges_mon,prob_mon