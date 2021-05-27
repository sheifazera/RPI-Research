# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 10:26:08 2021

@author: sheifazera
"""
import itertools
import random
import numpy as np
import time
from prettytable import PrettyTable, ALL
import networkx as nx
from pyomo.environ import *
from shortestpath_networkx import create_asymmetric_uncertainty_shortest_path_interdiction_nx
import matplotlib.pyplot as plt  


random.seed(a=631996)


def create_connected_graph(N,E):
    G=nx.random_tree(N)
    if nx.is_connected(G)==False:
        raise Exception ('Not a spanning tree!')
    for (i,j) in G.edges: #Adds data to spanning tree portion of graph
        p=random.random()
        q=0.5 *p
        G.edges[i,j]['length'] = -np.log(p)
        G.edges[i,j]['interdicted_length']=-np.log(q)+np.log(p)
    G_comp=nx.complete_graph(N)
    while G.number_of_edges() < 0.5* E: #generates remaining edges with data
        edge=random.sample(list(set(G_comp.edges)-set(G.edges)),1)
        for (i,j) in edge:
            p=random.random()
            q=0.5*p
            G.add_edge(i,j,length=-np.log(p), interdicted_length=-np.log(q)+np.log(p))
    G=G.to_directed()
    if nx.is_strongly_connected(G)==False:
        raise Exception('Directed graph is not strongly connected')
    nx.draw(G, with_labels=True)
    plt.show()
    plt.savefig("path.png")
    return G
 
def create_R(G):
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
    return R
    


time_10=[]
time_100=[]
time_1000=[]    
B_10=3
B_100=30
B_1000=300 
vdim_10=30
vdim_100=300
vdim_1000=3000
opt=SolverFactory('gurobi')
opt.options['TimeLimit']=1200
    
for loop in range(1,21):
    print(f'Starting loop {loop}')
    table=PrettyTable()
    table.hrules=ALL
    table.field_names=["Number of Nodes", "Model Build Time", "Solve Time"]
    
    
    
    G_10=create_connected_graph(10,30)
    R_10=create_R(G_10)
    print('Finished creating the graph')
    s=random.choice(list(G_10.nodes))
    t=random.choice(list(set(G_10.nodes)-{s}))
    start=time.time()
    M_10=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G_10,R_10,s,t,vdim_10,B_10)
    end=time.time()
    model_creation_time=end-start
    print(f'Finished creating the model in {model_creation_time} seconds')
    start=time.time()
    opt.solve(M_10)
    end=time.time()
    solve_time=end-start
    print(f'Finished solving 10 node model in {solve_time} seconds')
    table.add_row([10,model_creation_time,solve_time])
    time_10.append(solve_time)
    
    
    
    G_100=create_connected_graph(100,300)
    R_100=create_R(G_100)
    print('Finished creating the graph')
    s=random.choice(list(G_100.nodes))
    t=random.choice(list(set(G_100.nodes)-{s}))
    start=time.time()
    M_100=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G_100,R_100,s,t,vdim_100,B_100)
    end=time.time()
    model_creation_time=end-start
    print(f'Finished creating the model in {model_creation_time} seconds')
    start=time.time()
    opt.solve(M_100)
    end=time.time()
    solve_time=end-start
    print(f'Finished solving 10 node model in {solve_time} seconds')
    table.add_row([100,model_creation_time,solve_time])
    time_100.append(solve_time)
    
    
    
    
    G_1000=create_connected_graph(1000,3000)
    R_1000=create_R(G_1000)
    print('Finished creating the graph')
    s=random.choice(list(G_1000.nodes))
    t=random.choice(list(set(G_1000.nodes)-{s}))
    start=time.time()
    M_1000=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G_1000,R_1000,s,t,vdim_1000,B_1000)
    end=time.time()
    model_creation_time=end-start
    print(f'Finished creating the model in {model_creation_time} seconds')
    start=time.time()
    opt.solve(M_1000)
    end=time.time()
    solve_time=end-start
    print(f'Finished solving 10 node model in {solve_time} seconds')
    table.add_row([1000,model_creation_time,solve_time])
    time_1000.append(solve_time)
    
    
    table_txt=table.get_string()
    with open('TimeTrials_tables.txt','a') as file:
        file.write(table_txt)
    

    
    
    
    
    
        
        




