# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 13:13:37 2021

@author: spunlag
"""
import networkx as nx
import time
import matplotlib.pyplot as plt
from shortestpath_networkx import create_asymmetric_uncertainty_shortest_path_interdiction_nx, create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx, return_paths,  return_paths_multiple_ST
import random
import numpy as np
from pyomo.environ import *



D={(0,1):(1,0.5), (0,2):(0.8,0.4), (1,3):(1,1),(2,3):(1,1)}
G=nx.DiGraph()

for (i,j) in D.keys():
    (p,q)=D[(i,j)]
    G.add_edge(i,j)
    G[i][j]['length']=-np.log(p)
    G[i][j]['interdicted_length']=-np.log(q)+np.log(p)
    
vdim=len(G.edges)
R={}
    
for (i,j) in set(G.edges):
    r=G[i][j]['interdicted_length']
    k=0
    for (u,v) in set(G.edges):
        if i==u and j==v:
            R[(i,j,k)]=r
        else:
            R[(i,j,k)]=0
        k=k+1
      
        
S=0
T=3
B=2
M=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R,S,T,vdim,B)

opt=SolverFactory('gurobi_direct')
opt.solve(M)
M.x.pprint()
M.w.pprint()
M.z.pprint()

#%% Regret Avoided= Improvement, Regret=Objective Price You Pay

def get_IR(M,G,R,s,t,vdim,B):
    #Solve the nominal model
    R0={}
    vdim0=0
    M_nom=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R0, s,t,vdim0,B)
    opt.solve(M_nom)
    #(paths,length)=return_paths(M_nom,G,R0,s,t)
    #print('Nominal SPI: Path')
    #print(paths)
    #print(length)
    #Solve a model that has the fixed interdictions from the nominal model but adds in robustness
    M_nom_rob=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R, s,t,vdim,B)
    M_nom_rob.mods=ConstraintList()
    #print('Nominal SPI: Interdicted Arcs')
    for (i,j) in G.edges:
        if value(M_nom.x[(i,j)]) >= 0.9:
            M_nom_rob.mods.add(M_nom_rob.x[(i,j)]==1)
            #print(f'({i},{j})')
        else:
            M_nom_rob.mods.add(M_nom_rob.x[(i,j)]==0)
    opt.solve(M_nom_rob)
    #M_nom_rob.x.pprint()
    #(paths,length)=return_paths(M_nom_rob,G,R,s,t)
    #print('Robust Path Given Nominal Interdictions')
    #print(paths)
    #print(length)
    #Make comparisons
    #print(value(M.Obj))
    #print(value(M_nom.Obj))
    #print(value(M_nom_rob.Obj))
    Regret=100*(exp(-value(M.Obj))-exp(-value(M_nom.Obj)))/exp(-value(M.Obj))
    Improvement=100*(exp(-value(M_nom_rob.Obj))-exp(-value(M.Obj)))/exp(-value(M_nom_rob.Obj))

    return (Improvement,Regret)

(RegretAvoided,R)=get_IR(M,G,R,S,T,vdim,B)
print(RegretAvoided)