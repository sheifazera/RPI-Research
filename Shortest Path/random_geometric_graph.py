# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 09:19:13 2021

@author: sheifazera
"""

# -*- coding: utf-8 -*-

"""

Created on Thu May 27 09:33:14 2021

 

@author: spunlag

"""

 

import numpy as np

import itertools as it

import matplotlib.pyplot as plt

import random

import time

from pyomo.environ import *

#from prettytable import PrettyTable, ALL

from shortestpath_networkx import create_asymmetric_uncertainty_shortest_path_interdiction_nx, return_paths

import networkx as nx

import pickle

 

opt=SolverFactory('gurobi_direct')

opt.options['TimeLimit']=1800

 

Budgets={5,10,15,20}

Nodes={50,100,150,200}

random.seed(a=631996)

#%%

data={}

for loop in range(1,6):

    for N in Nodes:

        connected=False

        while connected is False:

            G=nx.random_geometric_graph(N, 1.5/sqrt(N), dim=2, p=2)

            connected=nx.is_connected(G)

        G=G.to_directed()

       

        for (i,j) in G.edges:

           p=random.random()

           q=random.uniform(0.25*p, 0.75*p)

           G[i][j]['length'] = -np.log(p)

           G[i][j]['interdicted_length']=-np.log(q)+np.log(p)

       

        R={}

        vdim=len(G.edges)

       

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

       

        for B in Budgets:

            M=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R,s,t,vdim,B)

            start=time.time()

            opt.solve(M)

            end=time.time()

            data[(N,B,loop)]=(end-start)

           

        print(f'Finished with {N} nodes')

    print(f'Finished with loop {loop}')

 
filename='data_1_5'
outfile = open(filename,'wb')
pickle.dump(data,outfile)
outfile.close()
     

#%%

pickle_off = open("data_1_5","rb")
data = pickle.load(pickle_off)

for loop in range(6,11):

   for N in Nodes:

       connected=False

       while connected is False:

           G=nx.random_geometric_graph(N, 1.5/sqrt(N), dim=2, p=2)

           connected=nx.is_connected(G)

       G=G.to_directed()

      

       for (i,j) in G.edges:

          p=random.random()

          q=random.uniform(0.25*p, 0.75*p)

          G[i][j]['length'] = -np.log(p)

          G[i][j]['interdicted_length']=-np.log(q)+np.log(p)

      

       R={}

       vdim=len(G.edges)

      

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

      

       for B in Budgets:

           M=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R,s,t,vdim,B)

           start=time.time()

           opt.solve(M)

           end=time.time()

           data[(N,B,loop)]=(end-start)

          

       print(f'Finished with {N} nodes')

   print(f'Finished with loop {loop}')  

    

filename='data_1_10'
outfile = open(filename,'wb')
pickle.dump(data,outfile)
outfile.close()  

 

#%%
pickle_off = open("data_1_10","rb")
data = pickle.load(pickle_off)

for loop in range(11,16):

    for N in Nodes:

        connected=False

        while connected is False:

            G=nx.random_geometric_graph(N, 1.5/sqrt(N), dim=2, p=2)

            connected=nx.is_connected(G)

        G=G.to_directed()

       

        for (i,j) in G.edges:

           p=random.random()

           q=random.uniform(0.25*p, 0.75*p)

           G[i][j]['length'] = -np.log(p)

           G[i][j]['interdicted_length']=-np.log(q)+np.log(p)

       

        R={}

        vdim=len(G.edges)

       

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

       

        for B in Budgets:

            M=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R,s,t,vdim,B)

            start=time.time()

            opt.solve(M)

            end=time.time()

            data[(N,B,loop)]=(end-start)

           

        print(f'Finished with {N} nodes')

    print(f'Finished with loop {loop}')

         

filename='data_1_15'
outfile = open(filename,'wb')
pickle.dump(data,outfile)
outfile.close()         

#%%

pickle_off = open("data_1_15","rb")
data = pickle.load(pickle_off) 

for loop in range(16,21):

    for N in Nodes:

        connected=False

        while connected is False:

            G=nx.random_geometric_graph(N, 1.5/sqrt(N), dim=2, p=2)

            connected=nx.is_connected(G)

        G=G.to_directed()

       

        for (i,j) in G.edges:

           p=random.random()

           q=random.uniform(0.25*p, 0.75*p)

           G[i][j]['length'] = -np.log(p)

           G[i][j]['interdicted_length']=-np.log(q)+np.log(p)

       

        R={}

        vdim=len(G.edges)

       

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

       

        for B in Budgets:

            M=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R,s,t,vdim,B)

            start=time.time()

            opt.solve(M)

            end=time.time()

            data[(N,B,loop)]=(end-start)

           

        print(f'Finished with {N} nodes')

    print(f'Finished with loop {loop}')  

   

filename='data_1_20'
outfile = open(filename,'wb')
pickle.dump(data,outfile)
outfile.close()      

#%% Data Visualization
pickle_off = open("data_1_20","rb")
data = pickle.load(pickle_off) 
 

for N in Nodes:

    fig1, ax1 = plt.subplots()

    plt.title('%i Nodes' % (N))

    #ax1.set_title('Solve time')

    data_array=np.empty([20,4])

    for loop in range(1,21):

        for B in Budgets:

            data_array[loop-1,int(B/5 -1)]=data[(N,B,loop)]

    ax1.boxplot(data_array, labels=[5,10,15,20], showfliers=False)

    ax1.set_xlabel('Budget')

    ax1.set_ylabel('Solve Time (s)')   

            

            

for B in Budgets:

    fig1, ax1 = plt.subplots()

    plt.title('Budget = %i' % (B))          

    data_array=np.empty([20,4])

    for loop in range(1,21):

        for N in Nodes:

            data_array[loop-1,int(N/50 -1)]=data[(N,B,loop)]

    ax1.boxplot(data_array, labels=[50,100,150,200], showfliers=False)

    ax1.set_xlabel('Number of Nodes')

    ax1.set_ylabel('Solve Time (s)')          