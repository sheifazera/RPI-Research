# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 10:19:46 2021

@author: spunlag
"""

import networkx as nx
import time
import matplotlib.pyplot as plt
from shortestpath_networkx import create_asymmetric_uncertainty_shortest_path_interdiction_nx
import random
import numpy as np
from pyomo.environ import *
import pickle


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
#%%
matrix = np.loadtxt('siouxfalls.txt', usecols=(0,1,3))

Prob = {(int(a),int(b)) : c for a,b,c in matrix}
for (i,j) in Prob.keys():
    if Prob[(i,j)] ==10:
        Prob[(i,j)] = (0.1,0.05) #p_ij, q_ij 
    elif Prob[(i,j)] ==9 :
        Prob[(i,j)] = (0.2, 0.1)
    elif Prob[(i,j)] ==8 :  
        Prob[(i,j)] = (0.3,0.15)
    elif Prob[(i,j)] ==7 :  
        Prob[(i,j)] = (0.4,0.2)
    elif Prob[(i,j)] ==6 :  
        Prob[(i,j)] = (0.5, 0.25)
    elif Prob[(i,j)] ==5 :  
        Prob[(i,j)] = (0.6, 0.3)
    elif Prob[(i,j)] ==4 :  
        Prob[(i,j)] = (0.7, 0.35)
    elif Prob[(i,j)] ==3 :  
        Prob[(i,j)] = (0.8, 0.4)
    elif Prob[(i,j)] ==2 :  
        Prob[(i,j)] = (0.9, 0.45)    
    else:
        raise Exception('Original Sioux Falls value out of range 2-10')

random.seed(a=631996)

G=nx.DiGraph()

for (i,j) in Prob.keys():
    (p,_)=Prob[(i,j)]
    q=random.uniform(0.25*p, 0.75*p)
    G.add_edge(i,j)
    G[i][j]['length']=-np.log(p)
    G[i][j]['interdicted_length']=-np.log(q)+np.log(p)
    
vdim=len(G.edges)

Budgets={2,4,6,8}
r_fracs=list(np.linspace(0.25,1,4))
s=24
t=10

for r_frac in r_fracs:
    R={}
    for (i,j) in set(G.edges):
        r=G[i][j]['interdicted_length']
        k=0
        for (u,v) in set(G.edges):
            if i==u and j==v:
                R[(i,j,k)]=r_frac*r
            else:
                R[(i,j,k)]=0
            k=k+1
    
    for B in Budgets: 
        M=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R,s,t,vdim,B)
        opt=SolverFactory('gurobi_direct')
        opt.solve(M)
        (Regret_Avoided, Regret)=get_IR(M,G,R,s,t,vdim,B)
        print(f'r_frac={r_frac} and B={B}')
        print(f'Regret Avoided={Regret_Avoided}')
    