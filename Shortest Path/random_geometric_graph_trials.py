# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 14:11:53 2021

@author: sheifazera
"""

import numpy as np
import itertools as it
import matplotlib.pyplot as plt
import random
import time
from pyomo.environ import *
from prettytable import PrettyTable, ALL
from shortestpath_networkx import create_asymmetric_uncertainty_shortest_path_interdiction_nx, return_paths
import networkx as nx
import pickle


#%%
random.seed(a=631996)
opt=SolverFactory('gurobi')
opt.options['TimeLimit']=1800

data={}
N=100
for loop in range(1,21):
    connected=False
    while connected is False:
        G=nx.random_geometric_graph(N,1.5/sqrt(N), dim=2, p=2)
        connected=nx.is_connected(G)
    plt.subplot(121)
    nx.draw(G, with_labels=False, font_weight='bold', node_size=25)
    G=G.to_directed()
    
    for (i,j) in G.edges:
       p=random.random()
       q=0.5 *p
       G.edges[i,j]['length'] = -np.log(p)
       G.edges[i,j]['interdicted_length']=-np.log(q)+np.log(p)
    
    s=random.choice(list(G.nodes))
    t=random.choice(list(set(G.nodes)-{s}))
    
    paths=nx.shortest_path(G,s,t,'length')
    length=len(paths)
    vdim=len(list(G.edges))
    B=length
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
    print('Finished creating model')
    start=time.time()
    opt.solve(M)
    end=time.time()
    (paths, lengths)=return_paths(M,G,R,s,t)
    number_of_paths=len(paths)
    data[(N,loop)]=(end-start, number_of_paths, paths[0], lengths, len(paths[0]))
    print(f'Finished loop {loop}')
   #%% 
N=500
for loop in range(1,21):
    connected=False
    while connected is False:
        G=nx.random_geometric_graph(N,1.5/sqrt(N), dim=2, p=2)
        connected=nx.is_connected(G)
    G=G.to_directed()
    for (i,j) in G.edges:
       p=random.random()
       q=0.5 *p
       G.edges[i,j]['length'] = -np.log(p)
       G.edges[i,j]['interdicted_length']=-np.log(q)+np.log(p)
    
    
    s=random.choice(list(G.nodes))
    t=random.choice(list(set(G.nodes)-{s}))
    
    paths=nx.shortest_path(G,s,t,'length')
    length=len(paths)
    
    vdim=len(list(G.edges))
    B=length
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
    start=time.time()
    opt.solve(M)
    end=time.time()
    (paths, lengths)=return_paths(M,G,R,s,t)
    number_of_paths=len(paths)
    data[(N,loop)]=(end-start, number_of_paths, paths[0], lengths, len(paths[0]))   
    print(f'Finished loop {loop}')
N=1000
for loop in range(1,21):
    connected=False
    while connected is False:
        G=nx.random_geometric_graph(N,1.5/sqrt(N), dim=2, p=2)
        connected=nx.is_connected(G)
    G=G.to_directed()
    for (i,j) in G.edges:
       p=random.random()
       q=0.5 *p
       G.edges[i,j]['length'] = -np.log(p)
       G.edges[i,j]['interdicted_length']=-np.log(q)+np.log(p)
    
    s=random.choice(list(G.nodes))
    t=random.choice(list(set(G.nodes)-{s}))
    
    paths=nx.shortest_path(G,s,t,'length')
    length=len(paths)
    
    vdim=len(list(G.edges))
    B=length
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
    start=time.time()
    opt.solve(M)
    end=time.time()
    (paths, lengths)=return_paths(M,G,R,s,t)
    number_of_paths=len(paths)
    data[(N,loop)]=(end-start, number_of_paths, paths[0], lengths, len(paths[0]))
    print(f'Finished loop {loop}')
#%%   
N=5000
for loop in range(1,21):
    connected=False
    while connected is False:
        G=nx.random_geometric_graph(N,1.5/sqrt(N), dim=2, p=2)
        connected=nx.is_connected(G)
    G=G.to_directed()
    for (i,j) in G.edges:
       p=random.random()
       q=0.5 *p
       G.edges[i,j]['length'] = -np.log(p)
       G.edges[i,j]['interdicted_length']=-np.log(q)+np.log(p)
    
    
    vdim=len(list(G.edges))
    B=int(.25*vdim)
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
    
    s=random.choice(list(G.nodes))
    t=random.choice(list(set(G.nodes)-{s}))
    
    M=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R,s,t,vdim,B)
    start=time.time()
    opt.solve(M)
    end=time.time()
    (paths, lengths)=return_paths(M,G,R,s,t)
    number_of_paths=len(paths)
    data[(N,loop)]=(end-start, number_of_paths, paths[0], lengths[0], len(paths[0]))
    
N=10000
for loop in range(1,21):
    connected=False
    while connected is False:
        G=nx.random_geometric_graph(N,1.5/sqrt(N), dim=2, p=2)
        connected=nx.is_connected(G)
    G=G.to_directed()
    for (i,j) in G.edges:
       p=random.random()
       q=0.5 *p
       G.edges[i,j]['length'] = -np.log(p)
       G.edges[i,j]['interdicted_length']=-np.log(q)+np.log(p)
    
    
    vdim=len(list(G.edges))
    B=int(.25*vdim)
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
    
    s=random.choice(list(G.nodes))
    t=random.choice(list(set(G.nodes)-{s}))
    
    M=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R,s,t,vdim,B)
    start=time.time()
    opt.solve(M)
    end=time.time()
    (paths, lengths)=return_paths(M,G,R,s,t)
    number_of_paths=len(paths)
    data[(N,loop)]=(end-start, number_of_paths, paths[0], lengths[0], len(paths[0]))


#%% Data Visualization
data_100_time=np.empty([20,1])
N=100
for i in range(1,21):
    (times, number_of_paths, path, length, number_of_nodes)=data[(N,i)]
    data_100_time[i-1,0]=times

data_500_time=np.empty([20,1])
N=500
for i in range(1,21):
    (times, number_of_paths, path, length, number_of_nodes)=data[(N,i)]
    data_100_time[i-1,0]=times

data_1000_time=np.empty([20,1])
N=1000
for i in range(1,21):
    (times, number_of_paths, path, length, number_of_nodes)=data[(N,i)]
    data_100_time[i-1,0]=times


collect=[data_100_time,data_500_time, data_1000_time]


fig1, ax1 = plt.subplots()
ax1.set_title('Solve time')
ax1.boxplot(collect, labels=[100, 500,1000], showfliers=False)
ax1.set_xlabel('Number of Nodes')
ax1.set_ylabel('Solve Time (s)')

fig1, ax1 = plt.subplots()
ax1.set_title('Solve time')
ax1.boxplot(collect, labels=[100,500,1000], showfliers=True)
ax1.set_xlabel('Number of Nodes')
ax1.set_ylabel('Solve Time (s)')
