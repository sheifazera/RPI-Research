# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 09:25:12 2021

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
from shortestpath_networkx import return_paths_multiple_ST , create_asymmetric_uncertainty_shortest_path_interdiction_nx, create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx #, return_path_robust
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

I=set(it.combinations(list(G.nodes),2))


Budgets={5,10,15,20,25,30,35,40}

#%%
data_time_array=np.empty([8,4,10])
data_prob_e_array=np.empty([8,4,10])
data_prob_d_array=np.empty([8,4,10])
#data{(budget, size of instance, trial):(time to solve, probability of evasion, probability of evasion to evader)}
random.seed(a=631996)
trials=10
opt=SolverFactory('gurobi')
opt.options['TimeLimit']=1200
for i in range(1,trials+1):
    I_16=set(random.sample(I, 16))
    I_8=set(random.sample(I_16, 8))
    I_4=set(random.sample(I_8, 4))
    I_2=set(random.sample(I_4, 2))
    
    for B in Budgets:
        
        
        M=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_16,vdim,B)
        start=time.time()
        opt.solve(M)
        end=time.time()
        SPD=100000
        SPE=0
        for (S,T) in I_16:
            (paths,length)=return_paths_multiple_ST(M,G,R,S,T)
            if length <SPD:
                SPD=length
                SPE=M.STBlocks[(S,T)].d[T].value
            
        data_time_array[int(B/5-1)][3][i-1]=end-start
        data_prob_e_array[int(B/5-1)][3][i-1]=SPD
        data_prob_d_array[int(B/5-1)][3][i-1]=SPE
        
        #I8
        
        
        M=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
        start=time.time()
        opt.solve(M)
        end=time.time()
        SPD=100000
        SPE=0
        for (S,T) in I_8:
            (paths,length)=return_paths_multiple_ST(M,G,R,S,T)
            if length <SPD:
                SPD=length
                SPE=M.STBlocks[(S,T)].d[T].value
            
        data_time_array[int(B/5-1)][2][i-1]=end-start
        data_prob_e_array[int(B/5-1)][2][i-1]=SPD
        data_prob_d_array[int(B/5-1)][2][i-1]=SPE
        
        #I4
        
        
        M=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_4,vdim,B)
        start=time.time()
        opt.solve(M)
        end=time.time()
        SPD=100000
        SPE=0
        for (S,T) in I_4:
            (paths,length)=return_paths_multiple_ST(M,G,R,S,T)
            if length <SPD:
                SPD=length
                SPE=M.STBlocks[(S,T)].d[T].value
            
        data_time_array[int(B/5-1)][1][i-1]=end-start
        data_prob_e_array[int(B/5-1)][1][i-1]=SPD
        data_prob_d_array[int(B/5-1)][1][i-1]=SPE
        #I2
        
        
        
        M=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_2,vdim,B)
        start=time.time()
        opt.solve(M)
        end=time.time()
        SPD=100000
        SPE=0
        for (S,T) in I_2:
            (paths,length)=return_paths_multiple_ST(M,G,R,S,T)
            if length <SPD:
                SPD=length
                SPE=M.STBlocks[(S,T)].d[T].value
            
        data_time_array[int(B/5-1)][0][i-1]=end-start
        data_prob_e_array[int(B/5-1)][0][i-1]=SPD
        data_prob_d_array[int(B/5-1)][0][i-1]=SPE
        print(f' For budget {B}, I have finished loop {i}')



filename='Budget_MST_time'
outfile = open(filename,'wb') 
pickle.dump(data_time_array,outfile)
outfile.close()     

filename='Budget_MST_prob_e'
outfile = open(filename,'wb') 
pickle.dump(data_prob_e_array,outfile)
outfile.close()     

filename='Budget_MST_prob_d'
outfile = open(filename,'wb') 
pickle.dump(data_prob_d_array,outfile)
outfile.close()        
        
#%% Data Averaging and Line plot
pickle_off = open("Budget_MST_time","rb")
data_time_array=pickle.load(pickle_off)

pickle_off = open("Budget_MST_prob_e","rb")
data_prob_e_array=pickle.load(pickle_off)

pickle_off = open("Budget_MST_prob_d","rb")
data_prob_d_array=pickle.load(pickle_off)


data_avg_e=np.empty([8,4])
data_avg_d=np.empty([8,4])

B_label=np.array([5,10,15,20,25,30,35,40])
for k in range(0,10):
    plt.figure()
    plt.plot(B_label, np.exp(-data_prob_e_array[:,0,k]), label = "2 (S,T) Pairs")
    plt.plot(B_label, np.exp(-data_prob_e_array[:,1,k]), label = "4 (S,T) Pairs")
    plt.plot(B_label, np.exp(-data_prob_e_array[:, 2,k]), label = "8 (S,T) Pairs")
    plt.plot(B_label, np.exp(-data_prob_e_array[:, 3,k]), label = "16 (S,T) Pairs")
    plt.xlabel('Budget')
    plt.ylabel('Probability of Evasion')
    plt.title('Budget, Number of (S,T) Pairs and Defender Evasion')
    plt.legend() 
    
    plt.figure()
    plt.plot(B_label, np.exp(-data_prob_d_array[:,0,k]), label = "2 (S,T) Pairs")
    plt.plot(B_label, np.exp(-data_prob_d_array[:,1,k]), label = "4 (S,T) Pairs")
    plt.plot(B_label, np.exp(-data_prob_d_array[:, 2,k]), label = "8 (S,T) Pairs")
    plt.plot(B_label, np.exp(-data_prob_d_array[:, 3,k]), label = "16 (S,T) Pairs")
    plt.xlabel('Budget')
    plt.ylabel('Probability of Evasion')
    plt.title('Budget, Number of (S,T) Pairs and Evader Evasion')
    plt.legend() 
    
    plt.figure()
    plt.plot(B_label, np.exp(-data_time_array[:,0,k]), label = "2 (S,T) Pairs")
    plt.plot(B_label, np.exp(-data_time_array[:,1,k]), label = "4 (S,T) Pairs")
    plt.plot(B_label, np.exp(-data_time_array[:, 2,k]), label = "8 (S,T) Pairs")
    plt.plot(B_label, np.exp(-data_time_array[:, 3,k]), label = "16 (S,T) Pairs")
    plt.xlabel('Budget')
    plt.ylabel('Time to Solve')
    plt.title('Budget, Number of (S,T) Pairs and Solve Time')
    plt.legend() 
    