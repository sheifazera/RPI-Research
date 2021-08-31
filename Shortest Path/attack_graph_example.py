# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 08:15:43 2021

@author: spunlag
"""

import networkx as nx
import time
import matplotlib.pyplot as plt
from shortestpath_networkx import create_asymmetric_uncertainty_shortest_path_interdiction_nx, create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx, return_paths,  return_paths_multiple_ST
import random
import numpy as np
from pyomo.environ import *



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
    print(exp(-value(M_nom_rob.Obj)))
    Regret=100*(exp(-value(M.Obj))-exp(-value(M_nom.Obj)))/exp(-value(M.Obj))
    Improvement=100*(exp(-value(M_nom_rob.Obj))-exp(-value(M.Obj)))/exp(-value(M_nom_rob.Obj))

    return (Improvement,Regret)


#%% Dictionary Definition

D={(0,101): (.99,.417), (0,102):(0.883,0.25), (101,1): (1,1), (102,1): (1,1), (1,121):(0.99,0.917), (1,122):(0.99, 0.5), (121,2): (1,1), (122,2): (1,1), (2,3):(0.917, 0.883), (3,341): (0.99, 0.883), (3,342): (0.99,0.917), (341,4): (1,1), (342,4): (1,1), (3,351): (0.99,0.883), (3,352): (0.99, 0.917), (351,5):(1,1), (352,5):(1,1), (4,461): (0.99, 0.75), (4,462): (0.99,.583), (461,6):(1,1), (462,6):(1,1), (4,471): (0.99, 0.75), (4,472): (0.99,.583), (471,7):(1,1), (472,7):(1,1), (4,481): (0.99, 0.75), (4,482): (0.99,.583), (481,8):(1,1), (482,8):(1,1), (4,491): (0.99, 0.75), (4,492): (0.99,.583), (491,9):(1,1), (492,9):(1,1), (4,4101): (0.99, 0.75), (4,4102): (0.99,.583), (4101,10):(1,1), (4102,10):(1,1), (5,5111): (0.99, 0.75), (5,5112): (0.99,.583), (5111,11):(1,1), (5112,11):(1,1),(5,5121): (0.99, 0.75), (5,5122): (0.99,.583), (5121,12):(1,1), (5122,12):(1,1), (5,5131): (0.99, 0.75), (5,5132): (0.99,.583), (5131,13):(1,1), (5132,13):(1,1),(5,5141): (0.99, 0.75), (5,5142): (0.99,.583), (5141,14):(1,1), (5142,14):(1,1)} 

D_SS=D
D_SS[(6,100)]=(1,1)
D_SS[(7,100)]=(1,1)
D_SS[(8,100)]=(1,1)
D_SS[(9,100)]=(1,1)
D_SS[(10,100)]=(1,1)
D_SS[(11,100)]=(1,1)
D_SS[(12,100)]=(1,1)
D_SS[(13,100)]=(1,1)
D_SS[(14,100)]=(1,1) 
#%% Super Sink Version

G=nx.DiGraph()

for (i,j) in D_SS.keys():
    (p,q)=D_SS[(i,j)]
    G.add_edge(i,j)
    G[i][j]['length']=-np.log(p)
    G[i][j]['interdicted_length']=-np.log(q)+np.log(p)

#nx.draw(G, with_labels=True, font_weight='bold')


R={}
vdim=0
B=12

'''
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
'''
s=0
t=100
M=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R,s,t,vdim,B)
load=125
M.min_load_shed=Constraint( expr = 125*(M.z[(9,100)] + M.w[(9,100)]) + 100*(M.z[(13,100)] + M.w[(13,100)]) + 90*(M.z[(12,100)] + M.w[(12,100)]) >= load )
opt=SolverFactory('gurobi_direct')
opt.solve(M)
print('Interdicted Arcs')
for (i,j) in G.edges():
    if value(M.x[(i,j)])>=.9:
        print(f'({i},{j})')
(paths, lengths)=return_paths(M,G,R,s,t)
print(paths)
print(exp(-lengths))
print(f'Nominal={exp(-M.d[t].value)}')
(Improvement,Regret)=get_IR(M,G,R,s,t,vdim,B)
print(Improvement)
print(Regret)
#%% Multiple (S,T) Version
G=nx.DiGraph()

for (i,j) in D.keys():
    (p,q)=D[(i,j)]
    G.add_edge(i,j)
    G[i][j]['length']=-np.log(p)
    G[i][j]['interdicted_length']=-np.log(q)+np.log(p)

#nx.draw(G, with_labels=True, font_weight='bold')


'''
R={}
vdim=0

'''
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
    


B=12
s=0
I={(s,6), (s,7), (s,8), (s,9), (s,10), (s,11), (s,12),(s,13), (s,14)}
M=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I,vdim,B)
opt=SolverFactory('gurobi_direct')
opt.solve(M)
print('Interdicted Arcs')
for (i,j) in G.edges():
    if value(M.x[(i,j)])>=.9:
        print(f'({i},{j})')
for (s,t) in I:     
    (paths,lengths)= return_paths_multiple_ST(M,G,R,s,t)
    print(f'({s},{t}):')
    print(paths)
    print(exp(-lengths))
    print(f'{exp(-M.STBlocks[(s,t)].d[t].value)}')
    (Improvement,Regret)=get_IR(M,G,R,s,t,vdim,B)
    print(Improvement)
    print(Regret)
    

#%% Symmetry Breaking Probability Modification


Budgets={2,4,6,8,10,12,14}
c1_list=np.linspace(0.01*0.833,0.51*0.833,6)
c2_list=np.linspace(0.01*0.25,0.51*0.25,6)
c_list=list(zip(c1_list, c2_list))
D={(0,101): (.99,.417), (0,102):(0.883,0.25), (101,1): (1,1), (102,1): (1,1), (1,121):(0.99,0.917), (1,122):(0.99, 0.5), (121,2): (1,1), (122,2): (1,1), (2,3):(0.917, 0.883), (3,341): (0.99, 0.883), (3,342): (0.99,0.917), (341,4): (1,1), (342,4): (1,1), (3,351): (0.99,0.883), (3,352): (0.99, 0.917), (351,5):(1,1), (352,5):(1,1),  (461,6):(1,1), (462,6):(1,1),  (471,7):(1,1), (472,7):(1,1),  (481,8):(1,1), (482,8):(1,1),  (491,9):(1,1), (492,9):(1,1),  (4101,10):(1,1), (4102,10):(1,1),  (5111,11):(1,1), (5112,11):(1,1), (5121,12):(1,1), (5122,12):(1,1),  (5131,13):(1,1), (5132,13):(1,1), (5141,14):(1,1), (5142,14):(1,1)} 
D_PLC={(4,461): (0.99, 0.75, 0), (4,462): (0.99,.583,0),  (4,471): (0.99, 0.75,0), (4,472): (0.99,.583,0),  (4,481): (0.99, 0.75,0), (4,482): (0.99,.583,0), (4,491): (0.99, 0.75,125), (4,492): (0.99,.583,125), (4,4101): (0.99, 0.75,0), (4,4102): (0.99,.583,0), (5,5111): (0.99, 0.75,0), (5,5112): (0.99,.583,0), (5,5121): (0.99, 0.75,90), (5,5122): (0.99,.583,90), (5,5131): (0.99, 0.75,100), (5,5132): (0.99,.583,100), (5,5141): (0.99, 0.75,0), (5,5142): (0.99,.583,0)}

D_SS=D
D_SS[(6,100)]=(1,1)
D_SS[(7,100)]=(1,1)
D_SS[(8,100)]=(1,1)
D_SS[(9,100)]=(1,1)
D_SS[(10,100)]=(1,1)
D_SS[(11,100)]=(1,1)
D_SS[(12,100)]=(1,1)
D_SS[(13,100)]=(1,1)
D_SS[(14,100)]=(1,1) 


G=nx.DiGraph()

for (i,j) in D_SS.keys():
    (p,q)=D_SS[(i,j)]
    G.add_edge(i,j)
    G[i][j]['length']=-np.log(p)
    G[i][j]['interdicted_length']=-np.log(q)+np.log(p)

for (c1,c2) in c_list:

    for (i,j) in D_PLC.keys():
        (p,q,LS)=D_PLC[(i,j)]
        G.add_edge(i,j)
        G[i][j]['length']=-np.log(p+c1*(LS/125))
        G[i][j]['interdicted_length']=-np.log(q+c2*(LS/125))+np.log(p+c1*(LS/125))
    
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
    #nx.draw(G, with_labels=True, font_weight='bold')
    for B in Budgets:
        s=0
        t=100
        M=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R,s,t,vdim,B)
        opt=SolverFactory('gurobi_direct')
        opt.solve(M)
        print(f'Budget={B}')
        print(f'c1={c1} and c2={c2}')
        print('Interdicted Arcs:')
        for (i,j) in set(G.edges):
            if value(M.x[(i,j)])>=0.9:
                print(f'({i},{j})')
        (paths, lengths)=return_paths(M,G,R,s,t)
        print(paths)
        print(exp(-lengths))

#%%
Budgets={2,4,6,8,10,12,14}
Budgets={2}
c1=0.01*0.99
c2=0.01*0.583
s=0
t=100


D={(0,101): (.99,.417), (0,102):(0.883,0.25), (101,1): (1,1), (102,1): (1,1), (1,121):(0.99,0.917), (1,122):(0.99, 0.5), (121,2): (1,1), (122,2): (1,1), (2,3):(0.917, 0.883), (3,341): (0.99, 0.883), (3,342): (0.99,0.917), (341,4): (1,1), (342,4): (1,1), (3,351): (0.99,0.883), (3,352): (0.99, 0.917), (351,5):(1,1), (352,5):(1,1),  (461,6):(1,1), (462,6):(1,1),  (471,7):(1,1), (472,7):(1,1),  (481,8):(1,1), (482,8):(1,1),  (491,9):(1,1), (492,9):(1,1),  (4101,10):(1,1), (4102,10):(1,1),  (5111,11):(1,1), (5112,11):(1,1), (5121,12):(1,1), (5122,12):(1,1),  (5131,13):(1,1), (5132,13):(1,1), (5141,14):(1,1), (5142,14):(1,1)} 
D_PLC={(4,461): (0.99, 0.75, 0), (4,462): (0.99,.583,0),  (4,471): (0.99, 0.75,0), (4,472): (0.99,.583,0),  (4,481): (0.99, 0.75,0), (4,482): (0.99,.583,0), (4,491): (0.99, 0.75,125), (4,492): (0.99,.583,125), (4,4101): (0.99, 0.75,0), (4,4102): (0.99,.583,0), (5,5111): (0.99, 0.75,0), (5,5112): (0.99,.583,0), (5,5121): (0.99, 0.75,90), (5,5122): (0.99,.583,90), (5,5131): (0.99, 0.75,100), (5,5132): (0.99,.583,100), (5,5141): (0.99, 0.75,0), (5,5142): (0.99,.583,0)}

D_SS=D
D_SS[(6,100)]=(1,1)
D_SS[(7,100)]=(1,1)
D_SS[(8,100)]=(1,1)
D_SS[(9,100)]=(1,1)
D_SS[(10,100)]=(1,1)
D_SS[(11,100)]=(1,1)
D_SS[(12,100)]=(1,1)
D_SS[(13,100)]=(1,1)
D_SS[(14,100)]=(1,1) 


G=nx.DiGraph()

for (i,j) in D_SS.keys():
    (p,q)=D_SS[(i,j)]
    G.add_edge(i,j)
    G[i][j]['length']=-np.log(p)
    G[i][j]['interdicted_length']=-np.log(q)+np.log(p)

G_nom=G

for (i,j) in D_PLC.keys():
        (p,q,LS)=D_PLC[(i,j)]
        G.add_edge(i,j)
        G_nom.add_edge(i,j)
        G_nom[i][j]['length']=-np.log(p)
        G[i][j]['length']=-np.log(p+c1*(LS/125))
        G_nom[i][j]['interdicted_length']=-np.log(q)+np.log(p)
        G[i][j]['interdicted_length']=-np.log(q+c2*(LS/125))+np.log(p+c1*(LS/125))
        
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

#vdim=0
#R={}

for B in Budgets:
    M=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R,s,t,vdim,B)
    opt=SolverFactory('gurobi_direct')
    opt.solve(M)
    print(f'B={B}')
    print('Interdicted Arcs')
    for (i,j) in G.edges():
        if value(M.x[(i,j)])>=.9:
            print(f'({i},{j})')
    (paths, lengths)=return_paths(M,G_nom,R,s,t)
    print(paths)
    print(exp(-lengths))
    print(f'Nominal={exp(-M.d[t].value)}')
    (Improvement,Regret)=get_IR(M,G,R,s,t,vdim,B)
    print(f'Regret Avoided={Improvement}%')


