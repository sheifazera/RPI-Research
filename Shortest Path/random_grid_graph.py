# -*- coding: utf-8 -*-

import networkx as nx
import time
import matplotlib.pyplot as plt
from shortestpath_networkx import create_asymmetric_uncertainty_shortest_path_interdiction_nx
import random
import numpy as np
from pyomo.environ import *
import pickle
#solverpathexe= "C:/Users\spunlag\Anaconda3\pkgs\gurobi-9.1.2-py38_0\gurobi_cl.exe"
#opt=SolverFactory('gurobi', executable=solverpathexe)
opt=SolverFactory('gurobi_direct')
opt.options['TimeLimit']=1800

random.seed(a=631996)
#There are m*n nodes in the graph
#(r,c) connects to (r+1,c), (r,c+1), (r+1,c+1)
 
def generate_grid_graph(n,m):
    G=nx.DiGraph()
    for r in range(1,n+1):
        for c in range(1,m+1):
        #For node number (r-1)+c
            
            if r!=n and c!=m:
                #Add an edge to node (r,c+1)
                G.add_edge((r-1)*m+c,(r-1)*m+c+1)
                #Add an edge to node (r+1,c)
                G.add_edge((r-1)*m+c,r*m+c)
                #Add an edge to node (r+1,c+1)
                G.add_edge((r-1)*m+c, r*m+c+1)
            elif r==n and c!=m:
                G.add_edge((r-1)*m+c,(r-1)*m+c+1)
            elif r!=n and c==m:
                G.add_edge((r-1)*m+c, r*m+c)
            #if r==n and c ==m, we have reached the end
    return G

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

Budgets={5,10,15,20}
Networks={(6,8), (6,10), (6,12), (8,8), (8,10), (8,12), (10,10)}

data={} #{(n,m,B):(time,%I, %R) %I: (prob of evasion for this model)-(prob of evasion for the nominal model decisions with robustness added after))/this model
#%R: (obj of this model)-(obj of nominal model) / (obj of this model)
for (n,m) in Networks:
    s=1
    t=n*m
    G=generate_grid_graph(n,m)
    for (i,j) in G.edges():
        p=random.random()
        q=random.uniform(0.25*p, 0.75*p)
        G[i][j]['length']=-np.log(p)
        G[i][j]['interdicted_length']=-np.log(q) + np.log(p)
        
        
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
    
    for B in Budgets:
        M=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R,s,t,vdim,B)
        start=time.time()
        results=opt.solve(M)
        end=time.time()
        if (results.solver.status == SolverStatus.ok) and (results.solver.termination_condition == TerminationCondition.optimal):
            (Improvement,Regret)=get_IR(M,G,R,s,t,vdim,B)
            data[(n,m,B)]=(end-start, Improvement,Regret)
        else:
            data[(n,m,B)]=(1800,0,0)
    print(f'Finished with Network ({n},{m})')
    
    
#%%
filename='Grid_graph'
outfile = open(filename,'wb')
pickle.dump(data,outfile)
outfile.close() 
 #%% Data Summary/Table Generation/Visualizations?

#save_game_load = 'Grid_graph' + '.pickle'
pickle_in = open("C:/Users/spunlag/Documents/Code/Grid_graph","rb")    


#%% 
random.seed(631996)
n=2
m=3
G=generate_grid_graph(2,3)
s=1
t=n*m
B=2
for (i,j) in G.edges():
        p=random.random()
        q=random.uniform(0.25*p, 0.75*p)
        G[i][j]['length']=-np.log(p)
        G[i][j]['interdicted_length']=-np.log(q) + np.log(p)
        
nx.draw(G, with_labels=True, font_weight='bold')        
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

M=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R,s,t,vdim,B)
tic=time.time()
opt.solve(M)
tok=time.time()
print('Asymmetric Information SPI')
(paths, length)=return_paths(M,G,R,s,t)
print(paths)
for (i,j) in G.edges():
    if value(M.x[(i,j)]) >=0.9:
        print(f'({i},{j})')
(Improvement,Regret)=get_IR(M,G,R,s,t,vdim,B)




 
        
        




