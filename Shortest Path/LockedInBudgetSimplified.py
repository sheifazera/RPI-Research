# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 09:03:56 2021

@author: sheifazera
"""

import numpy as np
import itertools as it
import random
import time
from pyomo.environ import *
from prettytable import PrettyTable, ALL
import networkx as nx
import pickle
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from shortestpath_networkx import return_interdicted_arcs, return_paths_multiple_ST , create_asymmetric_uncertainty_shortest_path_interdiction_nx, create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx #, return_path_robust
#%%
random.seed(a=631996)
matrix = np.loadtxt('eastern_mass.txt', usecols=(0,1,3))
Prob = {(int(a),int(b)) : c for a,b,c in matrix}

A=max(Prob.values())
B=min(Prob.values())
a=1
b=0.1

for (i,j) in Prob.keys():
    Prob[(i,j)]=(a-b)*(Prob[(i,j)]-B)/(A-B) + b
    
G=nx.DiGraph()
for (i,j) in Prob.keys():
    p=Prob[(i,j)]
    #q=0.5*p
    q=random.uniform(0.25*p, 0.75*p)
    G.add_edge(i,j,length=-np.log(p),interdicted_length=-np.log(q)+np.log(p))
    #r=0.25*q
    #G.add_edge(i,j,length=-np.log(p),interdicted_length=-np.log(q)+np.log(p), robust=-np.log(r)+np.log(q))

vdim=len(set(G.edges()))
R={}

for (i,j) in set(G.edges):
    #r=G[i][j]['robust']
    r=G[i][j]['interdicted_length']
    k=0
    for (u,v) in set(G.edges):
        if i==u and j==v:
            R[(i,j,k)]=r
        else:
            R[(i,j,k)]=0
        k=k+1    
#matrix=None
#Prob=None  

#First attempt
I_8={(4,36), (55,15), (74,22), (68, 2), (56,37), (33,49), (44,49), (70,5)}
I_8={(4,36), (55,15), (74,22), (68, 2), (56,37), (33,49), (44,49), (70,5), (1,15), (10,43), (55,32), (67,24), (53,21), (12,37), (5,40), (73,22)}
#Second attempt
#I_8={(4,36), (55,15), (68,2), (56,37), (33,49), (70,5), (51, 20), (62, 73) }
Budgets=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
#Budgets=[2,4,6,8,10,12,14,16,18,20]
#Budgets=[5,10,15,20]
bsize=len(Budgets)

#Experiment
#R={}
#vdim=0
#%% Function Definitions

def modify_optimistic(M_L, M_S, G):
    M_S.mods=ConstraintList()
    for (i,j) in set(G.edges):
        if M_L.x[(i,j)].value <= 0.9:
            M_S.mods.add(M_S.x[(i,j)]==0)
    return M_S

def modify_pessimistic(M_L, M_S, G):
    M_L.mods=ConstraintList()
    for (i,j) in set(G.edges):
        if M_S.x[(i,j)].value >= 0.9:
            M_L.mods.add(M_L.x[(i,j)]==1)
    return M_L

def collect_data(M,G,R,I):
    SPD=100000
    SPE=0
    for (S,T) in I:
        (paths,length)=return_paths_multiple_ST(M,G,R,S,T)
        if length <SPD:
            SPD=length
            SPE=M.STBlocks[(S,T)].d[T].value
            path=paths

    return (exp(-SPE),exp(-SPD),path)

opt=SolverFactory('gurobi')
opt.options['TimeLimit']=1200

#%% Pessimistic
data_ref=np.empty([bsize,1])
data_pess=np.empty([bsize,1])
i=0
for B in Budgets:
    print(f'Budget={B}')
    
    M_ref=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
    opt.solve(M_ref)
    i_list=return_interdicted_arcs(M_ref,G)
    print('Reference:')
    print(i_list)
    #get probability of evasion
    (d_e,d_d, path)=collect_data(M_ref,G,R,I_8)
    data_ref[(i,0)]=d_d
    print(path)
    
    
    if B>=Budgets[1]:
        #do modified version
        M_pess=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
        M_pess=modify_pessimistic(M_pess, M_pess_old, G)
        opt.solve(M_pess)
        i_list=return_interdicted_arcs(M_pess,G)
        (d_e,d_d, path)=collect_data(M_pess,G,R,I_8)
        data_pess[(i,0)]=d_d
        print('Pessimisitc:')
        print(i_list)
        print(path)
        M_pess_old=M_pess
        
        
    else:
        data_pess[(i,0)]=data_ref[(i,0)]
        M_pess_old=M_ref

    i=i+1
    
    
#%% Optimistic

data_op=np.empty([bsize,1])
Budgets.reverse()
i=0
for B in Budgets:
    print(f'Budget={B}')

    if B<=Budgets[1]:
        #do modified version
        M_op=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
        M_op=modify_optimistic(M_op_old, M_op, G)
        opt.solve(M_op)
        i_list=return_interdicted_arcs(M_op,G)
        (d_e,d_d, path)=collect_data(M_op,G,R,I_8)
        data_op[(i,0)]=d_d
        print(i_list)
        print(path)
        M_op_old=M_op
        
    else:
        M_ref=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
        opt.solve(M_ref)
        i_list=return_interdicted_arcs(M_ref,G)
        print(i_list)
        #get probability of evasion
        (d_e,d_d, path)=collect_data(M_ref,G,R,I_8)
        print(path)
        data_op[(i,0)]=d_d
        M_op_old=M_ref
    i=i+1   
    
Budgets.reverse()  
data_op=np.flip(data_op)
#%% Data Visualization
B_label=np.asarray(Budgets)
plt.figure()
plt.xlabel('Budget')
plt.ylabel('Robust Probability of Evasion')
plt.title('16 (S,T) Pairs')  
plt.plot(B_label, data_pess, label="Pessimistic", linestyle="-.", lw=4)
plt.plot(B_label, data_ref, label= "Reference", linestyle="--", marker="o", lw=4) 
plt.plot(B_label, data_op, label= "Optimistic", linestyle=":", lw=4)
 
plt.xticks(B_label)
plt.legend()