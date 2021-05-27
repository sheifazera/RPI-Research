# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 09:25:32 2021

@author: sheifazera
"""


import numpy as np
import itertools as it
import matplotlib.pyplot as plt
import random
import time
from pyomo.environ import *
from prettytable import PrettyTable, ALL
from shortestpath_networkx import create_shortest_path_nx, create_shortest_path_interdiction_nx, create_asymmetric_uncertainty_shortest_path_interdiction_nx, return_paths, similar_length, create_shortest_path_interdiction_nx
import networkx as nx
import pickle
import matplotlib.cm as cm

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
#%%
I=set(it.combinations(list(G.nodes),2))

pickle_off = open("ST_data_with_similar","rb")
data = pickle.load(pickle_off)


random.seed(a=631996)
bad_ST=[]
for (S,T) in data.keys():
    (_,_,_,t,_,_,_,_)=data[(S,T)]
    if t > 1:
        bad_ST.append((S,T))

#Sample len(bad_ST) not-bad pairs as well
good_ST=random.sample(list(I-set(bad_ST)),len(bad_ST))

ST=good_ST + bad_ST
#%% 
data_small={} #data[(S,T)]: {solve time, length of path, nodes in path, number of similar paths}
cutoff=1.1
Budgets={5, 7, 10, 12, 15, 17, 20}
loop=0
opt=SolverFactory('gurobi')
opt.options['TimeLimit']=1200
for B in Budgets:
    for (S,T) in ST:
        M=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R,S,T,vdim,B)   
        start=time.time()
        opt.solve(M)
        end=time.time()
        (paths, lengths)=return_paths(M,G,R,S,T)
        (similar_paths, no_of_similar_paths)=similar_length(M,G,R,S,T, cutoff)
        data_small[(S,T,B)]=(end-start, lengths, len(paths[0]), no_of_similar_paths)
        loop=loop+1
        print(f'{100*loop/(7*888)} percent complete')
    print(f'Finished with budget {B}')   
#%% 

    
filename='ST_data_subset'
outfile = open(filename,'wb')
pickle.dump(data_small,outfile)
outfile.close()


#%%
pickle_off = open("ST_data_subset","rb")
data_small = pickle.load(pickle_off)
#%% Data Sorting
data_5_bad=[]
data_5_good=[]
data_7_bad=[]
data_7_good=[]
data_10_bad=[]
data_10_good=[]
data_12_bad=[]
data_12_good=[]
data_15_bad=[]
data_15_good=[]
data_17_bad=[]
data_17_good=[]
data_20_bad=[]
data_20_good=[]

for (S,T,B) in data_small.keys():
    (t, l, n, s)=data_small[(S,T,B)]
    if B== 5:
        if (S,T) in bad_ST:
            data_5_bad.append([t,l,n,s])
        else:
            data_5_good.append([t,l,n,s])
    elif B==7:
        if (S,T) in bad_ST:
            data_7_bad.append([t,l,n,s])
        else:
            data_7_good.append([t,l,n,s])
    elif B==10:
        if (S,T) in bad_ST:
            data_10_bad.append([t,l,n,s])
        else:
            data_10_good.append([t,l,n,s])
    elif B==12:
        if (S,T) in bad_ST:
            data_12_bad.append([t,l,n,s])
        else:
            data_12_good.append([t,l,n,s])
    elif B==15:
        if (S,T) in bad_ST:
            data_15_bad.append([t,l,n,s])
        else:
            data_15_good.append([t,l,n,s])
    elif B==17:
        if (S,T) in bad_ST:
            data_17_bad.append([t,l,n,s])
        else:
            data_17_good.append([t,l,n,s])
    elif B==20:
        if (S,T) in bad_ST:
            data_20_bad.append([t,l,n,s])
        else:
            data_20_good.append([t,l,n,s])
        
    

data_5_bad_array=np.array(data_5_bad)
data_7_bad_array=np.array(data_7_bad)
data_10_bad_array=np.array(data_10_bad)
data_12_bad_array=np.array(data_12_bad)
data_15_bad_array=np.array(data_15_bad)
data_17_bad_array=np.array(data_17_bad)
data_20_bad_array=np.array(data_20_bad)

data_5_good_array=np.array(data_5_good)
data_7_good_array=np.array(data_7_good)
data_10_good_array=np.array(data_10_good)
data_12_good_array=np.array(data_12_good)
data_15_good_array=np.array(data_15_good)
data_17_good_array=np.array(data_17_good)
data_20_good_array=np.array(data_20_good)

#%% Get the B=25 data
data_25_good=[]
data_25_bad=[]
for (S,T) in good_ST:
    (_,_,_,t,_,l,n,s)=data[(S,T)]
    data_25_good.append([t,l,n,s])
for (S,T) in bad_ST:
    (_,_,_,t,_,l,n,s)=data[(S,T)]
    data_25_bad.append([t,l,n,s])
    
data_25_good_array=np.array(data_25_good)
data_25_bad_array=np.array(data_25_bad)


#%% Visualization

#Scatter Plots Length vs Time

fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(data_5_bad_array[:,1], data_5_bad_array[:,0],color='maroon')
ax.scatter(data_5_good_array[:,1], data_5_good_array[:,0], color='darkgreen')
ax.set_xlabel('Length of Path')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 5')


fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(data_7_bad_array[:,1], data_7_bad_array[:,0],color='firebrick')
ax.scatter(data_7_good_array[:,1], data_7_good_array[:,0], color='forestgreen')
ax.set_xlabel('Length of Path')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 7')


fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(data_10_bad_array[:,1], data_10_bad_array[:,0],color='indianred')
ax.scatter(data_10_good_array[:,1], data_10_good_array[:,0], color='olivedrab')
ax.set_xlabel('Length of Path')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 10')


fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(data_12_bad_array[:,1], data_12_bad_array[:,0],color='salmon')
ax.scatter(data_12_good_array[:,1], data_12_good_array[:,0], color='limegreen')
ax.set_xlabel('Length of Path')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 12')



fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(data_15_bad_array[:,1], data_15_bad_array[:,0],color='lightcoral')
ax.scatter(data_15_good_array[:,1], data_15_good_array[:,0], color='yellowgreen')
ax.set_xlabel('Length of Path')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 15')



fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(data_17_bad_array[:,1], data_17_bad_array[:,0],color='lightsalmon')
ax.scatter(data_17_good_array[:,1], data_17_good_array[:,0], color='lightgreen')
ax.set_xlabel('Length of Path')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 17')



fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(data_20_bad_array[:,1], data_20_bad_array[:,0],color='peachpuff')
ax.scatter(data_20_good_array[:,1], data_20_good_array[:,0], color='lightgreen')
ax.set_xlabel('Length of Path')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 20')


fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(data_25_bad_array[:,1], data_25_bad_array[:,0],color='mistyrose')
ax.scatter(data_25_good_array[:,1], data_25_good_array[:,0], color='palegreen')
ax.set_xlabel('Length of Path')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 25')

#%% 

#Scatter Plots Similar vs Time

fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(data_5_bad_array[:,3], data_5_bad_array[:,0],color='maroon')
ax.scatter(data_5_good_array[:,3], data_5_good_array[:,0], color='darkgreen')
ax.set_xlabel('Number of Similar Paths')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 5')


fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(data_7_bad_array[:,3], data_7_bad_array[:,0],color='firebrick')
ax.scatter(data_7_good_array[:,3], data_7_good_array[:,0], color='forestgreen')
ax.set_xlabel('Number of Similar Paths')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 7')


fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(data_10_bad_array[:,3], data_10_bad_array[:,0],color='indianred')
ax.scatter(data_10_good_array[:,3], data_10_good_array[:,0], color='olivedrab')
ax.set_xlabel('Number of Similar Paths')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 10')


fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(data_12_bad_array[:,3], data_12_bad_array[:,0],color='salmon')
ax.scatter(data_12_good_array[:,3], data_12_good_array[:,0], color='limegreen')
ax.set_xlabel('Number of Similar Paths')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 12')



fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(data_15_bad_array[:,3], data_15_bad_array[:,0],color='lightcoral')
ax.scatter(data_15_good_array[:,3], data_15_good_array[:,0], color='yellowgreen')
ax.set_xlabel('Number of Similar Paths')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 15')



fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(data_17_bad_array[:,3], data_17_bad_array[:,0],color='lightsalmon')
ax.scatter(data_17_good_array[:,3], data_17_good_array[:,0], color='lightgreen')
ax.set_xlabel('Number of Similar Paths')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 17')



fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(data_20_bad_array[:,3], data_20_bad_array[:,0],color='peachpuff')
ax.scatter(data_20_good_array[:,3], data_20_good_array[:,0], color='lightgreen')
ax.set_xlabel('Number of Similar Paths')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 20')


fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(data_25_bad_array[:,3], data_25_bad_array[:,0],color='mistyrose')
ax.scatter(data_25_good_array[:,3], data_25_good_array[:,0], color='palegreen')
ax.set_xlabel('Number of Similar Paths')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 25')

#%% Histograms


#%% Boxplots

data_5_array=np.vstack((data_5_good_array, data_5_bad_array))
data_7_array=np.vstack((data_7_good_array, data_7_bad_array))
data_10_array=np.vstack((data_10_good_array, data_10_bad_array))
data_12_array=np.vstack((data_12_good_array, data_12_bad_array))
data_15_array=np.vstack((data_15_good_array, data_15_bad_array))
data_17_array=np.vstack((data_17_good_array, data_17_bad_array))
data_20_array=np.vstack((data_20_good_array, data_20_bad_array))
data_25_array=np.vstack((data_25_good_array, data_25_bad_array))
#%%
fig1, ax1 = plt.subplots()
ax1.set_title('Solve time')
ax1.boxplot(data_5_array[:,0], labels=[5], showfliers=False)
ax1.set_xlabel('Budget')
ax1.set_ylabel('Solve Time (s)')

fig1, ax1 = plt.subplots()
ax1.set_title('Solve time')
ax1.boxplot(data_7_array[:,0], labels=[7], showfliers=False)
ax1.set_xlabel('Budget')
ax1.set_ylabel('Solve Time (s)')

fig1, ax1 = plt.subplots()
ax1.set_title('Solve time')
ax1.boxplot(data_10_array[:,0], labels=[10], showfliers=False)
ax1.set_xlabel('Budget')
ax1.set_ylabel('Solve Time (s)')

fig1, ax1 = plt.subplots()
ax1.set_title('Solve time')
ax1.boxplot(data_12_array[:,0], labels=[12], showfliers=False)
ax1.set_xlabel('Budget')
ax1.set_ylabel('Solve Time (s)')

fig1, ax1 = plt.subplots()
ax1.set_title('Solve time')
ax1.boxplot(data_15_array[:,0], labels=[15], showfliers=False)
ax1.set_xlabel('Budget')
ax1.set_ylabel('Solve Time (s)')

fig1, ax1 = plt.subplots()
ax1.set_title('Solve time')
ax1.boxplot(data_17_array[:,0], labels=[17], showfliers=False)
ax1.set_xlabel('Budget')
ax1.set_ylabel('Solve Time (s)')

fig1, ax1 = plt.subplots()
ax1.set_title('Solve time')
ax1.boxplot(data_20_array[:,0], labels=[20], showfliers=False)
ax1.set_xlabel('Budget')
ax1.set_ylabel('Solve Time (s)')

#%% Heat Maps 
from matplotlib.colors import LogNorm

x_min = np.min(data_5_array[:,3]) 
x_max = np.max(data_5_array[:,3]) 
  
y_min = np.min(data_5_array[:,0]) 
y_max = np.max(data_5_array[:,0]) 
  
x_bins = np.linspace(x_min, x_max, 31) 
y_bins = np.linspace(y_min, y_max, 101)

fig, ax = plt.subplots(figsize =(10, 7)) 
plt.hist2d(data_5_array[:,3], data_5_array[:,0], bins =[x_bins, y_bins], norm=LogNorm())
                
ax.set_xlabel('Number Similar Length Paths (10%)')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 5')
plt.colorbar()



x_min = np.min(data_7_array[:,3]) 
x_max = np.max(data_7_array[:,3]) 
  
y_min = np.min(data_7_array[:,0]) 
y_max = np.max(data_7_array[:,0]) 
  
x_bins = np.linspace(x_min, x_max, 31) 
y_bins = np.linspace(y_min, y_max, 101)

fig, ax = plt.subplots(figsize =(10, 7)) 
plt.hist2d(data_7_array[:,3], data_7_array[:,0], bins =[x_bins, y_bins], norm=LogNorm())
                
ax.set_xlabel('Number Similar Length Paths (10%)')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 7')
plt.colorbar()

x_min = np.min(data_10_array[:,3]) 
x_max = np.max(data_10_array[:,3]) 
  
y_min = np.min(data_10_array[:,0]) 
y_max = np.max(data_10_array[:,0]) 
  
x_bins = np.linspace(x_min, x_max, 31) 
y_bins = np.linspace(y_min, y_max, 101)

fig, ax = plt.subplots(figsize =(10, 7)) 
plt.hist2d(data_10_array[:,3], data_10_array[:,0], bins =[x_bins, y_bins], norm=LogNorm())
                
ax.set_xlabel('Number Similar Length Paths (10%)')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 10')
plt.colorbar()


x_min = np.min(data_12_array[:,3]) 
x_max = np.max(data_12_array[:,3]) 
  
y_min = np.min(data_12_array[:,0]) 
y_max = np.max(data_12_array[:,0]) 
  
x_bins = np.linspace(x_min, x_max, 31) 
y_bins = np.linspace(y_min, y_max, 101)

fig, ax = plt.subplots(figsize =(10, 7)) 
plt.hist2d(data_12_array[:,3], data_12_array[:,0], bins =[x_bins, y_bins], norm=LogNorm())
                
ax.set_xlabel('Number Similar Length Paths (10%)')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 12')
plt.colorbar()


x_min = np.min(data_15_array[:,3]) 
x_max = np.max(data_15_array[:,3]) 
  
y_min = np.min(data_15_array[:,0]) 
y_max = np.max(data_15_array[:,0]) 
  
x_bins = np.linspace(x_min, x_max, 31) 
y_bins = np.linspace(y_min, y_max, 101)

fig, ax = plt.subplots(figsize =(10, 7)) 
plt.hist2d(data_15_array[:,3], data_15_array[:,0], bins =[x_bins, y_bins], norm=LogNorm())
                
ax.set_xlabel('Number Similar Length Paths (10%)')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 15')
plt.colorbar()


x_min = np.min(data_17_array[:,3]) 
x_max = np.max(data_17_array[:,3]) 
  
y_min = np.min(data_17_array[:,0]) 
y_max = np.max(data_17_array[:,0]) 
  
x_bins = np.linspace(x_min, x_max, 31) 
y_bins = np.linspace(y_min, y_max, 101)

fig, ax = plt.subplots(figsize =(10, 7)) 
plt.hist2d(data_17_array[:,3], data_17_array[:,0], bins =[x_bins, y_bins], norm=LogNorm())
                
ax.set_xlabel('Number Similar Length Paths (10%)')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 17')
plt.colorbar()



x_min = np.min(data_20_array[:,3]) 
x_max = np.max(data_20_array[:,3]) 
  
y_min = np.min(data_20_array[:,0]) 
y_max = np.max(data_20_array[:,0]) 
  
x_bins = np.linspace(x_min, x_max, 31) 
y_bins = np.linspace(y_min, y_max, 101)

fig, ax = plt.subplots(figsize =(10, 7)) 
plt.hist2d(data_20_array[:,3], data_20_array[:,0], bins =[x_bins, y_bins], norm=LogNorm())
                
ax.set_xlabel('Number Similar Length Paths (10%)')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 20')
plt.colorbar()


x_min = np.min(data_25_array[:,3]) 
x_max = np.max(data_25_array[:,3]) 
  
y_min = np.min(data_25_array[:,0]) 
y_max = np.max(data_25_array[:,0]) 
  
x_bins = np.linspace(x_min, x_max, 31) 
y_bins = np.linspace(y_min, y_max, 101)

fig, ax = plt.subplots(figsize =(10, 7)) 
plt.hist2d(data_25_array[:,3], data_25_array[:,0], bins =[x_bins, y_bins], norm=LogNorm())
                
ax.set_xlabel('Number Similar Length Paths (10%)')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 25')
plt.colorbar()

#%% Heat Map Shortest Path LEngth


x_min = np.min(data_5_array[:,1]) 
x_max = np.max(data_5_array[:,1]) 
  
y_min = np.min(data_5_array[:,0]) 
y_max = np.max(data_5_array[:,0]) 
  
x_bins = np.linspace(x_min, x_max, 31) 
y_bins = np.linspace(y_min, y_max, 101)

fig, ax = plt.subplots(figsize =(10, 7)) 
plt.hist2d(data_5_array[:,1], data_5_array[:,0], bins =[x_bins, y_bins], norm=LogNorm())
                
ax.set_xlabel('Shortest Path Length')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 5')
plt.colorbar()



x_min = np.min(data_7_array[:,1]) 
x_max = np.max(data_7_array[:,1]) 
  
y_min = np.min(data_7_array[:,0]) 
y_max = np.max(data_7_array[:,0]) 
  
x_bins = np.linspace(x_min, x_max, 31) 
y_bins = np.linspace(y_min, y_max, 101)

fig, ax = plt.subplots(figsize =(10, 7)) 
plt.hist2d(data_7_array[:,1], data_7_array[:,0], bins =[x_bins, y_bins], norm=LogNorm())
                
ax.set_xlabel('Shortest Path Length')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 7')
plt.colorbar()

x_min = np.min(data_10_array[:,1]) 
x_max = np.max(data_10_array[:,1]) 
  
y_min = np.min(data_10_array[:,0]) 
y_max = np.max(data_10_array[:,0]) 
  
x_bins = np.linspace(x_min, x_max, 31) 
y_bins = np.linspace(y_min, y_max, 101)

fig, ax = plt.subplots(figsize =(10, 7)) 
plt.hist2d(data_10_array[:,1], data_10_array[:,0], bins =[x_bins, y_bins], norm=LogNorm())
                
ax.set_xlabel('Shortest Path Length')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 10')
plt.colorbar()


x_min = np.min(data_12_array[:,1]) 
x_max = np.max(data_12_array[:,1]) 
  
y_min = np.min(data_12_array[:,0]) 
y_max = np.max(data_12_array[:,0]) 
  
x_bins = np.linspace(x_min, x_max, 31) 
y_bins = np.linspace(y_min, y_max, 101)

fig, ax = plt.subplots(figsize =(10, 7)) 
plt.hist2d(data_12_array[:,1], data_12_array[:,0], bins =[x_bins, y_bins], norm=LogNorm())
                
ax.set_xlabel('Shortest Path Length')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 12')
plt.colorbar()


x_min = np.min(data_15_array[:,1]) 
x_max = np.max(data_15_array[:,1]) 
  
y_min = np.min(data_15_array[:,0]) 
y_max = np.max(data_15_array[:,0]) 
  
x_bins = np.linspace(x_min, x_max, 31) 
y_bins = np.linspace(y_min, y_max, 101)

fig, ax = plt.subplots(figsize =(10, 7)) 
plt.hist2d(data_15_array[:,1], data_15_array[:,0], bins =[x_bins, y_bins], norm=LogNorm())
                
ax.set_xlabel('Shortest Path Length')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 15')
plt.colorbar()


x_min = np.min(data_17_array[:,1]) 
x_max = np.max(data_17_array[:,1]) 
  
y_min = np.min(data_17_array[:,0]) 
y_max = np.max(data_17_array[:,0]) 
  
x_bins = np.linspace(x_min, x_max, 31) 
y_bins = np.linspace(y_min, y_max, 101)

fig, ax = plt.subplots(figsize =(10, 7)) 
plt.hist2d(data_17_array[:,1], data_17_array[:,0], bins =[x_bins, y_bins], norm=LogNorm())
                
ax.set_xlabel('Shortest Path Length')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 17')
plt.colorbar()



x_min = np.min(data_20_array[:,1]) 
x_max = np.max(data_20_array[:,1]) 
  
y_min = np.min(data_20_array[:,0]) 
y_max = np.max(data_20_array[:,0]) 
  
x_bins = np.linspace(x_min, x_max, 31) 
y_bins = np.linspace(y_min, y_max, 101)

fig, ax = plt.subplots(figsize =(10, 7)) 
plt.hist2d(data_20_array[:,1], data_20_array[:,0], bins =[x_bins, y_bins], norm=LogNorm())
                
ax.set_xlabel('Shortest Path Length')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 20')
plt.colorbar()


x_min = np.min(data_25_array[:,1]) 
x_max = np.max(data_25_array[:,1]) 
  
y_min = np.min(data_25_array[:,0]) 
y_max = np.max(data_25_array[:,0]) 
  
x_bins = np.linspace(x_min, x_max, 31) 
y_bins = np.linspace(y_min, y_max, 101)

fig, ax = plt.subplots(figsize =(10, 7)) 
plt.hist2d(data_25_array[:,1], data_25_array[:,0], bins =[x_bins, y_bins], norm=LogNorm())
                
ax.set_xlabel('Shortest Path Length')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Budget 25')
plt.colorbar()