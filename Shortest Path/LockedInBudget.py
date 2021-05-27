# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 09:12:10 2021

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
    q=0.5*p
    G.add_edge(i,j,length=-np.log(p),interdicted_length=-np.log(q)+np.log(p))

vdim=len(set(G.edges()))
R={}

for (i,j) in set(G.edges):
    l=G[i][j]['length']
    r=G[i][j]['interdicted_length']
    k=0
    for (u,v) in set(G.edges):
        if i==u and j==v:
            R[(i,j,k)]=r
        else:
            R[(i,j,k)]=0
        k=k+1    
matrix=None
Prob=None  


#Pick interesting options here
'''
I=set(it.combinations(list(G.nodes),2))
random.seed(a=631996)
I_16=set(random.sample(I-set(G.edges), 16))
I_8=set(random.sample(I_16, 8))
I_4=set(random.sample(I_8, 4))
'''
I_16={(4,36), (55,15), (74,22), (68, 2), (56,37), (33,49), (44,49), (70,5), (30, 14), (49, 24), (48,30), (58,6), (54,29), (70,21), (40,16), (52,24)}
I_8={(4,36), (55,15), (74,22), (68, 2), (56,37), (33,49), (44,49), (70,5)}
I_4={(4,36), (55,15), (74,22), (68,2)}
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

    return (exp(-SPE),exp(-SPD))

opt=SolverFactory('gurobi')
opt.options['TimeLimit']=1200


#%% Not locked-in Reference

data_e_16=[]
data_d_16=[]

data_e_8=[]
data_d_8=[]

data_e_4=[]
data_d_4=[]

print('Reference')
print('B=30')
B=30
M_30_16=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_16,vdim,B)
opt.solve(M_30_16)
(d_e,d_d)=collect_data(M_30_16,G,R,I_16)
data_e_16.append(d_e)
data_d_16.append(d_d)
i_list=return_interdicted_arcs(M_30_16,G)
print(i_list)


M_30_8=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
opt.solve(M_30_8)
(d_e,d_d)=collect_data(M_30_8,G,R,I_8)
data_e_8.append(d_e)
data_d_8.append(d_d)
i_list=return_interdicted_arcs(M_30_8,G)
print(i_list)

M_30_4=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_4,vdim,B)
opt.solve(M_30_4)
(d_e,d_d)=collect_data(M_30_4,G,R,I_4)
data_e_4.append(d_e)
data_d_4.append(d_d)
i_list=return_interdicted_arcs(M_30_4,G)
print(i_list)

print('B=25')
B=25
M_25_16=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_16,vdim,B)
opt.solve(M_25_16)
(d_e,d_d)=collect_data(M_25_16,G,R,I_16)
data_e_16.append(d_e)
data_d_16.append(d_d)
i_list=return_interdicted_arcs(M_25_16,G)
print(i_list)

M_25_8=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
opt.solve(M_25_8)
(d_e,d_d)=collect_data(M_25_8,G,R,I_8)
data_e_8.append(d_e)
data_d_8.append(d_d)
i_list=return_interdicted_arcs(M_25_8,G)
print(i_list)


M_25_4=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_4,vdim,B)
opt.solve(M_25_4)
(d_e,d_d)=collect_data(M_25_4,G,R,I_4)
data_e_4.append(d_e)
data_d_4.append(d_d)
i_list=return_interdicted_arcs(M_25_4,G)
print(i_list)

print('B=20')
B=20
M_20_16=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_16,vdim,B)
opt.solve(M_20_16)
(d_e,d_d)=collect_data(M_20_16,G,R,I_16)
data_e_16.append(d_e)
data_d_16.append(d_d)
i_list=return_interdicted_arcs(M_20_16,G)
print(i_list)


M_20_8=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
opt.solve(M_20_8)
(d_e,d_d)=collect_data(M_20_8,G,R,I_8)
data_e_8.append(d_e)
data_d_8.append(d_d)
i_list=return_interdicted_arcs(M_20_8,G)
print(i_list)

M_20_4=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_4,vdim,B)
opt.solve(M_20_4)
(d_e,d_d)=collect_data(M_20_4,G,R,I_4)
data_e_4.append(d_e)
data_d_4.append(d_d)
i_list=return_interdicted_arcs(M_20_4,G)
print(i_list)

print('B=15')
B=15
M_15_16=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_16,vdim,B)
opt.solve(M_15_16)
(d_e,d_d)=collect_data(M_15_16,G,R,I_16)
data_e_16.append(d_e)
data_d_16.append(d_d)
i_list=return_interdicted_arcs(M_15_16,G)
print(i_list)

M_15_8=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
opt.solve(M_15_8)
(d_e,d_d)=collect_data(M_15_8,G,R,I_8)
data_e_8.append(d_e)
data_d_8.append(d_d)
i_list=return_interdicted_arcs(M_15_8,G)
print(i_list)

M_15_4=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_4,vdim,B)
opt.solve(M_15_4)
(d_e,d_d)=collect_data(M_15_4,G,R,I_4)
data_e_4.append(d_e)
data_d_4.append(d_d)
i_list=return_interdicted_arcs(M_15_4,G)
print(i_list)

print('B=10')
B=10
M_10_16=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_16,vdim,B)
opt.solve(M_10_16)
(d_e,d_d)=collect_data(M_10_16,G,R,I_16)
data_e_16.append(d_e)
data_d_16.append(d_d)
i_list=return_interdicted_arcs(M_10_16,G)
print(i_list)

M_10_8=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
opt.solve(M_10_8)
(d_e,d_d)=collect_data(M_10_8,G,R,I_8)
data_e_8.append(d_e)
data_d_8.append(d_d)
i_list=return_interdicted_arcs(M_10_8,G)
print(i_list)

M_10_4=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_4,vdim,B)
opt.solve(M_10_4)
(d_e,d_d)=collect_data(M_10_4,G,R,I_4)
data_e_4.append(d_e)
data_d_4.append(d_d)
i_list=return_interdicted_arcs(M_10_4,G)
print(i_list)

print('B=5')
B=5
M_5_16=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_16,vdim,B)
opt.solve(M_5_16)
(d_e,d_d)=collect_data(M_5_16,G,R,I_16)
data_e_16.append(d_e)
data_d_16.append(d_d)
i_list=return_interdicted_arcs(M_5_16,G)
print(i_list)

M_5_8=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
opt.solve(M_5_8)
(d_e,d_d)=collect_data(M_5_8,G,R,I_8)
data_e_8.append(d_e)
data_d_8.append(d_d)
i_list=return_interdicted_arcs(M_5_8,G)
print(i_list)

M_5_4=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_4,vdim,B)
opt.solve(M_5_4)
(d_e,d_d)=collect_data(M_5_4,G,R,I_4)
data_e_4.append(d_e)
data_d_4.append(d_d)
i_list=return_interdicted_arcs(M_5_4,G)
print(i_list)
#%% Optimistic

data_o_e_16=[]
data_o_d_16=[]

data_o_e_8=[]
data_o_d_8=[]

data_o_e_4=[]
data_o_d_4=[]

print('Optimistic')
print('B=30')
B=30

M_O_30_16=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_16,vdim,B)
opt.solve(M_O_30_16)
(d_e,d_d)=collect_data(M_O_30_16,G,R,I_16)
data_o_e_16.append(d_e)
data_o_d_16.append(d_d)
i_list=return_interdicted_arcs(M_O_30_16,G)
print(i_list)


M_O_30_8=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
opt.solve(M_O_30_8)
(d_e,d_d)=collect_data(M_O_30_8,G,R,I_8)
data_o_e_8.append(d_e)
data_o_d_8.append(d_d)
i_list=return_interdicted_arcs(M_O_30_8,G)
print(i_list)


M_O_30_4=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_4,vdim,B)
opt.solve(M_O_30_4)
(d_e,d_d)=collect_data(M_O_30_4,G,R,I_4)
data_o_e_4.append(d_e)
data_o_d_4.append(d_d)
i_list=return_interdicted_arcs(M_O_30_4,G)
print(i_list)

print('B=25')
B=25

M_O_25_16=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_16,vdim,B)
M_O_25_16=modify_optimistic(M_O_30_16, M_O_25_16, G)
opt.solve(M_O_25_16)
(d_e,d_d)=collect_data(M_O_25_16,G,R,I_16)
data_o_e_16.append(d_e)
data_o_d_16.append(d_d)
i_list=return_interdicted_arcs(M_O_25_16,G)
print(i_list)



M_O_25_8=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
M_O_25_8=modify_optimistic(M_O_30_8, M_O_25_8, G)
opt.solve(M_O_25_8)
(d_e,d_d)=collect_data(M_O_25_8,G,R,I_8)
data_o_e_8.append(d_e)
data_o_d_8.append(d_d)
i_list=return_interdicted_arcs(M_O_25_8,G)
print(i_list)


M_O_25_4=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_4,vdim,B)
M_O_25_4=modify_optimistic(M_O_30_4, M_O_25_4, G)
opt.solve(M_O_25_4)
(d_e,d_d)=collect_data(M_O_25_4,G,R,I_4)
data_o_e_4.append(d_e)
data_o_d_4.append(d_d)
i_list=return_interdicted_arcs(M_O_25_4,G)
print(i_list)

print('B=20')
B=20

M_O_20_16=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_16,vdim,B)
M_O_20_16=modify_optimistic(M_O_25_16, M_O_20_16, G)
opt.solve(M_O_20_16)
(d_e,d_d)=collect_data(M_O_20_16,G,R,I_16)
data_o_e_16.append(d_e)
data_o_d_16.append(d_d)
i_list=return_interdicted_arcs(M_O_20_16,G)
print(i_list)


M_O_20_8=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
M_O_20_8=modify_optimistic(M_O_25_8, M_O_20_8, G)
opt.solve(M_O_20_8)
(d_e,d_d)=collect_data(M_O_20_8,G,R,I_8)
data_o_e_8.append(d_e)
data_o_d_8.append(d_d)
i_list=return_interdicted_arcs(M_O_20_8,G)
print(i_list)


M_O_20_4=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_4,vdim,B)
M_O_20_4=modify_optimistic(M_O_25_4, M_O_20_4, G)
opt.solve(M_O_20_4)
(d_e,d_d)=collect_data(M_O_20_4,G,R,I_4)
data_o_e_4.append(d_e)
data_o_d_4.append(d_d)
i_list=return_interdicted_arcs(M_O_20_4,G)
print(i_list)

print('B=15')
B=15

M_O_15_16=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_16,vdim,B)
M_O_15_16=modify_optimistic(M_O_20_16, M_O_15_16, G)
opt.solve(M_O_15_16)
(d_e,d_d)=collect_data(M_O_15_16,G,R,I_16)
data_o_e_16.append(d_e)
data_o_d_16.append(d_d)
i_list=return_interdicted_arcs(M_O_15_16,G)
print(i_list)


M_O_15_8=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
M_O_15_8=modify_optimistic(M_O_20_8, M_O_15_8, G)
opt.solve(M_O_15_8)
(d_e,d_d)=collect_data(M_O_15_8,G,R,I_8)
data_o_e_8.append(d_e)
data_o_d_8.append(d_d)
i_list=return_interdicted_arcs(M_O_15_8,G)
print(i_list)


M_O_15_4=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_4,vdim,B)
M_O_15_4=modify_optimistic(M_O_20_4, M_O_15_4, G)
opt.solve(M_O_15_4)
(d_e,d_d)=collect_data(M_O_15_4,G,R,I_4)
data_o_e_4.append(d_e)
data_o_d_4.append(d_d)
i_list=return_interdicted_arcs(M_O_15_4,G)
print(i_list)

print('B=10')
B=10

M_O_10_16=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_16,vdim,B)
M_O_10_16=modify_optimistic(M_O_15_16, M_O_10_16, G)
opt.solve(M_O_10_16)
(d_e,d_d)=collect_data(M_O_10_16,G,R,I_16)
data_o_e_16.append(d_e)
data_o_d_16.append(d_d)
i_list=return_interdicted_arcs(M_O_10_16,G)
print(i_list)


M_O_10_8=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
M_O_10_8=modify_optimistic(M_O_15_8, M_O_10_8, G)
opt.solve(M_O_10_8)
(d_e,d_d)=collect_data(M_O_10_8,G,R,I_8)
data_o_e_8.append(d_e)
data_o_d_8.append(d_d)
i_list=return_interdicted_arcs(M_O_10_8,G)
print(i_list)


M_O_10_4=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_4,vdim,B)
M_O_10_4=modify_optimistic(M_O_15_4, M_O_10_4, G)
opt.solve(M_O_10_4)
(d_e,d_d)=collect_data(M_O_10_4,G,R,I_4)
data_o_e_4.append(d_e)
data_o_d_4.append(d_d)
i_list=return_interdicted_arcs(M_O_10_4,G)
print(i_list)
print('B=5')
B=5

M_O_5_16=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_16,vdim,B)
M_O_5_16=modify_optimistic(M_O_10_16, M_O_5_16, G)
opt.solve(M_O_5_16)
(d_e,d_d)=collect_data(M_O_5_16,G,R,I_16)
data_o_e_16.append(d_e)
data_o_d_16.append(d_d)
i_list=return_interdicted_arcs(M_O_5_16,G)
print(i_list)

M_O_5_8=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
M_O_5_8=modify_optimistic(M_O_10_8, M_O_5_8, G)
opt.solve(M_O_5_8)
(d_e,d_d)=collect_data(M_O_5_8,G,R,I_8)
data_o_e_8.append(d_e)
data_o_d_8.append(d_d)
i_list=return_interdicted_arcs(M_O_5_8,G)
print(i_list)

M_O_5_4=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_4,vdim,B)
M_O_5_4=modify_optimistic(M_O_10_4, M_O_5_4, G)
opt.solve(M_O_5_4)
(d_e,d_d)=collect_data(M_O_5_4,G,R,I_4)
data_o_e_4.append(d_e)
data_o_d_4.append(d_d)
i_list=return_interdicted_arcs(M_O_5_4,G)
print(i_list)


#%% Pessimistic
data_p_e_16=[]
data_p_d_16=[]

data_p_e_8=[]
data_p_d_8=[]

data_p_e_4=[]
data_p_d_4=[]

print('Pessimistic')
print('B=5')
B=5

M_P_5_16=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_16,vdim,B)
opt.solve(M_P_5_16)
(d_e,d_d)=collect_data(M_P_5_16,G,R,I_16)
data_p_e_16.append(d_e)
data_p_d_16.append(d_d)
i_list=return_interdicted_arcs(M_P_5_16,G)
print(i_list)

M_P_5_8=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
opt.solve(M_P_5_8)
(d_e,d_d)=collect_data(M_P_5_8,G,R,I_8)
data_p_e_8.append(d_e)
data_p_d_8.append(d_d)
i_list=return_interdicted_arcs(M_P_5_8,G)
print(i_list)

M_P_5_4=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_4,vdim,B)
opt.solve(M_P_5_4)
(d_e,d_d)=collect_data(M_P_5_4,G,R,I_4)
data_p_e_4.append(d_e)
data_p_d_4.append(d_d)
i_list=return_interdicted_arcs(M_P_5_4,G)
print(i_list)
print('B=10')
B=10

M_P_10_16=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_16,vdim,B)
M_P_10_16=modify_pessimistic(M_P_10_16, M_P_5_16, G)
opt.solve(M_P_10_16)
(d_e,d_d)=collect_data(M_P_10_16,G,R,I_16)
data_p_e_16.append(d_e)
data_p_d_16.append(d_d)
i_list=return_interdicted_arcs(M_P_10_16,G)
print(i_list)

M_P_10_8=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
M_P_10_8=modify_pessimistic(M_P_10_8, M_P_5_8, G)
opt.solve(M_P_10_8)
(d_e,d_d)=collect_data(M_P_10_8,G,R,I_8)
data_p_e_8.append(d_e)
data_p_d_8.append(d_d)
i_list=return_interdicted_arcs(M_P_10_8,G)
print(i_list)

M_P_10_4=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_4,vdim,B)
M_P_10_4=modify_pessimistic(M_P_10_4, M_P_5_4, G)
opt.solve(M_P_10_4)
(d_e,d_d)=collect_data(M_P_10_4,G,R,I_4)
data_p_e_4.append(d_e)
data_p_d_4.append(d_d)
i_list=return_interdicted_arcs(M_P_10_4,G)
print(i_list)


print('B=15')
B=15

M_P_15_16=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_16,vdim,B)
M_P_15_16=modify_pessimistic(M_P_15_16, M_P_10_16, G)
opt.solve(M_P_15_16)
(d_e,d_d)=collect_data(M_P_15_16,G,R,I_16)
data_p_e_16.append(d_e)
data_p_d_16.append(d_d)
i_list=return_interdicted_arcs(M_P_15_16,G)
print(i_list)



M_P_15_8=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
M_P_15_8=modify_pessimistic(M_P_15_8, M_P_10_8, G)
opt.solve(M_P_15_8)
(d_e,d_d)=collect_data(M_P_15_8,G,R,I_8)
data_p_e_8.append(d_e)
data_p_d_8.append(d_d)
i_list=return_interdicted_arcs(M_P_15_8,G)
print(i_list)



M_P_15_4=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_4,vdim,B)
M_P_15_4=modify_pessimistic(M_P_15_4, M_P_10_4, G)
opt.solve(M_P_15_4)
(d_e,d_d)=collect_data(M_P_15_4,G,R,I_4)
data_p_e_4.append(d_e)
data_p_d_4.append(d_d)
i_list=return_interdicted_arcs(M_P_15_4,G)
print(i_list)


print('B=20')
B=20

M_P_20_16=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_16,vdim,B)
M_P_20_16=modify_pessimistic(M_P_20_16, M_P_15_16, G)
opt.solve(M_P_20_16)
(d_e,d_d)=collect_data(M_P_20_16,G,R,I_16)
data_p_e_16.append(d_e)
data_p_d_16.append(d_d)
i_list=return_interdicted_arcs(M_P_20_16,G)
print(i_list)



M_P_20_8=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
M_P_20_8=modify_pessimistic(M_P_20_8, M_P_15_8, G)
opt.solve(M_P_20_8)
(d_e,d_d)=collect_data(M_P_20_8,G,R,I_8)
data_p_e_8.append(d_e)
data_p_d_8.append(d_d)
i_list=return_interdicted_arcs(M_P_20_8,G)
print(i_list)

M_P_20_4=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_4,vdim,B)
M_P_20_4=modify_pessimistic(M_P_20_4, M_P_15_4, G)
opt.solve(M_P_20_4)
(d_e,d_d)=collect_data(M_P_20_4,G,R,I_4)
data_p_e_4.append(d_e)
data_p_d_4.append(d_d)
i_list=return_interdicted_arcs(M_P_20_4,G)
print(i_list)
print('B=25')
B=25

M_P_25_16=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_16,vdim,B)
M_P_25_16=modify_pessimistic(M_P_25_16, M_P_20_16, G)
opt.solve(M_P_25_16)
(d_e,d_d)=collect_data(M_P_25_16,G,R,I_16)
data_p_e_16.append(d_e)
data_p_d_16.append(d_d)
i_list=return_interdicted_arcs(M_P_25_16,G)
print(i_list)


M_P_25_8=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
M_P_25_8=modify_pessimistic(M_P_25_8, M_P_20_8, G)
opt.solve(M_P_25_8)
(d_e,d_d)=collect_data(M_P_25_8,G,R,I_8)
data_p_e_8.append(d_e)
data_p_d_8.append(d_d)
i_list=return_interdicted_arcs(M_P_25_8,G)
print(i_list)

M_P_25_4=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_4,vdim,B)
M_P_25_4=modify_pessimistic(M_P_25_4, M_P_20_4, G)
opt.solve(M_P_25_4)
(d_e,d_d)=collect_data(M_P_25_4,G,R,I_4)
data_p_e_4.append(d_e)
data_p_d_4.append(d_d)
i_list=return_interdicted_arcs(M_P_25_4,G)
print(i_list)

print('B=30')
B=30

M_P_30_16=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_16,vdim,B)
M_P_30_16=modify_pessimistic(M_P_30_16, M_P_25_16, G)
opt.solve(M_P_30_16)
(d_e,d_d)=collect_data(M_P_30_16,G,R,I_16)
data_p_e_16.append(d_e)
data_p_d_16.append(d_d)
i_list=return_interdicted_arcs(M_P_30_16,G)
print(i_list)

M_P_30_8=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
M_P_30_8=modify_pessimistic(M_P_30_8, M_P_25_8, G)
opt.solve(M_P_30_8)
(d_e,d_d)=collect_data(M_P_30_8,G,R,I_8)
data_p_e_8.append(d_e)
data_p_d_8.append(d_d)
i_list=return_interdicted_arcs(M_P_30_8,G)
print(i_list)

M_P_30_4=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_4,vdim,B)
M_P_30_4=modify_pessimistic(M_P_30_4, M_P_25_4, G)
opt.solve(M_P_30_4)
(d_e,d_d)=collect_data(M_P_30_4,G,R,I_4)
data_p_e_4.append(d_e)
data_p_d_4.append(d_d)
i_list=return_interdicted_arcs(M_P_30_4,G)
print(i_list)
#%%  Results
#for each size:
#Plot budget versus probability of evasion
#three lines: reference, optimistic evasion


data_e_16=np.asarray(data_e_16)
data_d_16=np.asarray(data_d_16)

data_e_8=np.asarray(data_e_8)
data_d_8=np.asarray(data_d_8)

data_e_4=np.asarray(data_e_4)
data_d_4=np.asarray(data_d_4)

data_o_e_16=np.asarray(data_o_e_16)
data_o_d_16=np.asarray(data_o_d_16)

data_o_e_8=np.asarray(data_o_e_8)
data_o_d_8=np.asarray(data_o_d_8)

data_o_e_4=np.asarray(data_o_e_4)
data_o_d_4=np.asarray(data_o_d_4)


data_p_e_16.reverse()
data_p_d_16.reverse()

data_p_e_16=np.asarray(data_p_e_16)
data_p_d_16=np.asarray(data_p_d_16)

data_p_e_8.reverse()
data_p_d_8.reverse()
data_p_e_8=np.asarray(data_p_e_8)
data_p_d_8=np.asarray(data_p_d_8)

data_p_e_4.reverse() 
data_p_d_4.reverse()                     
data_p_e_4=np.asarray(data_p_e_4)
data_p_d_4=np.asarray(data_p_d_4)


B_label=np.array([5,10,15,20,25,30])
plt.figure()
plt.plot(B_label, data_e_16, label = "Reference")
plt.plot(B_label, data_o_e_16, label = "Optimistic")
plt.plot(B_label, data_p_e_16, label = "Pessimistic")
plt.xlabel('Budget')
plt.ylabel('Nominal Probability of Evasion')
plt.title('16 (S,T) Pairs')
plt.legend()

plt.figure()
plt.plot(B_label, data_e_8, label = "Reference")
plt.plot(B_label, data_o_e_8, label = "Optimistic")
plt.plot(B_label, data_p_e_8, label = "Pessimistic")
plt.xlabel('Budget')
plt.ylabel('Nominal Probability of Evasion')
plt.title('8 (S,T) Pairs')
plt.legend()

plt.figure()
plt.plot(B_label, data_e_4, label = "Reference")
plt.plot(B_label, data_o_e_4, label = "Optimistic")
plt.plot(B_label, data_p_e_4, label = "Pessimistic")
plt.xlabel('Budget')
plt.ylabel('Nominal Probability of Evasion')
plt.title('4 (S,T) Pairs')
plt.legend()


plt.figure()
plt.plot(B_label, data_d_16, label = "Reference")
plt.plot(B_label, data_o_d_16, label = "Optimistic")
plt.plot(B_label, data_p_d_16, label = "Pessimistic")
plt.xlabel('Budget')
plt.ylabel('Robust Probability of Evasion')
plt.title('16 (S,T) Pairs')
plt.legend()

plt.figure()
plt.plot(B_label, data_d_8, label = "Reference")
plt.plot(B_label, data_o_d_8, label = "Optimistic")
plt.plot(B_label, data_p_d_8, label = "Pessimistic")
plt.xlabel('Budget')
plt.ylabel('Robust Probability of Evasion')
plt.title('8 (S,T) Pairs')
plt.legend()

plt.figure()
plt.plot(B_label, data_d_4, label = "Reference")
plt.plot(B_label, data_o_d_4, label = "Optimistic")
plt.plot(B_label, data_p_d_4, label = "Pessimistic")
plt.xlabel('Budget')
plt.ylabel('Robust Probability of Evasion')
plt.title('4 (S,T) Pairs')
plt.legend()

#%% 
#separate optimistic and pessimistic
#plot budget versus regret
#all sizes on same plot

plt.figure()
plt.plot(B_label, abs(data_d_16 - data_o_d_16), label="16 (S,T) Pairs")
plt.plot(B_label, abs(data_d_8 - data_o_d_8), label="8 (S,T) Pairs")
plt.plot(B_label, abs(data_d_4 - data_o_d_4), label="4 (S,T) Pairs")
plt.xlabel('Budget')
plt.ylabel('Regret')
plt.title('Optmistic')
plt.legend()

plt.figure()
plt.plot(B_label, abs(data_d_16 - data_p_d_16), label="16 (S,T) Pairs")
plt.plot(B_label, abs(data_d_8 - data_p_d_8), label="8 (S,T) Pairs")
plt.plot(B_label, abs(data_d_4 - data_p_d_4), label="4 (S,T) Pairs")
plt.xlabel('Budget')
plt.ylabel('Regret')
plt.title('Pessimistic')
plt.legend()

#%% 


plt.figure()
plt.plot(B_label, abs(data_d_16 - data_o_d_16)/data_d_16, label="16 (S,T) Pairs")
plt.plot(B_label, abs(data_d_8 - data_o_d_8)/data_d_8, label="8 (S,T) Pairs")
plt.plot(B_label, abs(data_d_4 - data_o_d_4)/data_d_4, label="4 (S,T) Pairs")
plt.xlabel('Budget')
plt.ylabel('Normalized Regret')
plt.title('Optmistic')
plt.legend()

plt.figure()
plt.plot(B_label, abs(data_d_16 - data_p_d_16)/data_d_16, label="16 (S,T) Pairs")
plt.plot(B_label, abs(data_d_8 - data_p_d_8)/data_d_8, label="8 (S,T) Pairs")
plt.plot(B_label, abs(data_d_4 - data_p_d_4)/data_d_4, label="4 (S,T) Pairs")
plt.xlabel('Budget')
plt.ylabel('Normalized Regret')
plt.title('Pessimistic')
plt.legend()