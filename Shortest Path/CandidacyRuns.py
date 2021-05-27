# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 10:11:54 2021

@author: sheifazera
"""
'''
Re-running and expanding past tests for my candidacy paper
'''

import numpy as np
from pyomo.environ import *
import networkx as nx
from shortestpath import create_shortest_path, create_shortest_path_robust_ellipsoid, return_path, return_path_robust #, return_path_robust
from shortestpath_networkx import create_asymmetric_uncertainty_shortest_path_interdiction_nx, return_interdicted_arcs, return_paths
#%%
matrix = np.loadtxt('siouxfalls.txt', usecols=(0,1,3))
#print(matrix) 
#%% 
#Sioux Falls Uncertainty Factors
D = {(int(a),int(b)) : c for a,b,c in matrix}
N=set(range(1,25))
s=3
t=20

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

D={}
for (i,j) in Prob.keys():
    (p,q)=Prob[(i,j)]
    D[(i,j)]=(-np.log(p))
        
udim=7
N_disrupt={8,10,11,15,16,20,22}
P={}
i=0
for n in N_disrupt: #for all nodes disrupted (uncertainty dimension)
    for (u,v) in D.keys(): #check all edges 
        if u==n or v==n: #edge enters the disrupted node 
            P[(u,v,i)]=-np.log(0.5)
        elif (u,n) in Prob.keys() or (v,n) in Prob.keys(): #u or v is adjacent to the disrupted node
            P[(u,v,i)]=-np.log(0.25) 
        else:    
            P[(u,v,i)]=0
    i=i+1 
ME=create_shortest_path_robust_ellipsoid(D,N,P,s,t,udim)
opt=SolverFactory('gurobi')
results_robust_ellipsoid=opt.solve(ME) 
path_robust=return_path_robust(ME,D,N,P,s,t,udim)

print(path_robust)
ME.u.pprint()
#%%    
#Sioux Falls Square P

D = {(int(a),int(b)) : c for a,b,c in matrix}
N=set(range(1,25))
s=3
t=20

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
        
D={}
for (i,j) in Prob.keys():
    (p,q)=Prob[(i,j)]
    D[(i,j)]=(-np.log(p))

udim=len(D.keys())
P={}

for (i,j) in D.keys():
    k=0
    for (u,v) in D.keys():
        if i==u and j==v:
            P[(i,j,k)]=-np.log(0.25)
        else:
            P[(i,j,k)]=0
        k=k+1    
        
ME=create_shortest_path_robust_ellipsoid(D,N,P,s,t,udim)

results_robust_ellipsoid=opt.solve(ME) 
path_robust=return_path_robust(ME,D,N,P,s,t,udim)

print(path_robust)     
k=1   
for (i,j) in D.keys():
    if ME.u[k].value > 0.1:
        print(i,j)
        print(ME.u[k].value)
    k=k+1
#%%
P={}
udim=0
ME=create_shortest_path_robust_ellipsoid(D,N,P,s,t,udim)
results_robust_ellipsoid=opt.solve(ME) 
path_robust=return_path_robust(ME,D,N,P,s,t,udim)

print(path_robust)        
ME.u.pprint()        



#%% Asymmetric Uncertainty - Single (S,T) Pair
D = {(int(a),int(b)) : c for a,b,c in matrix}
N=set(range(1,25))
s=3
t=20

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
G=nx.DiGraph()
for (i,j) in Prob.keys():
    (p,q)=Prob[(i,j)]
    G.add_edge(i,j,length=-np.log(p),interdicted_length=-np.log(q)+np.log(p))
    #r=0.25*q
    #G.add_edge(i,j,length=-np.log(p),interdicted_length=-np.log(q)+np.log(p), robust=-np.log(r)+np.log(q))
vdim=len(set(G.edges()))
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

vdim=0
R={}

B=5        
M_ref=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R,s,t,vdim,B)     
opt.solve(M_ref)   
arcs=return_interdicted_arcs(M_ref,G)
(path,length)=return_paths(M_ref,G,R,s,t)
print(path)
print(exp(-length))
print(arcs)