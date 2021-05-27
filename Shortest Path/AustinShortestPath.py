# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 12:43:05 2021

@author: sheifazera
"""

import numpy as np
import random
import time
from pyomo.environ import *
from prettytable import PrettyTable, ALL
from shortestpath import create_shortest_path, create_shortest_path_robust_ellipsoid, return_path , create_asymmetric_uncertainty_shortest_path_interdiction, create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST #, return_path_robust
import matplotlib.pyplot as plt  


matrix = np.loadtxt('Barcelona.txt', usecols=(0,1,3))
#print(matrix) 
print('Matrix fully loaded')
Prob = {(int(a),int(b)) : c for a,b,c in matrix}
matrix=None
print('Matrix converted to dictionary')

N=set()
for (i,j) in Prob.keys():
    if i not in N:
        N.add(i)
    if j not in N:
        N.add(j)
#N=set(range(1,1021)) #check this

A=max(Prob.values())
B=min(Prob.values())
a=1
b=0.1

for (i,j) in Prob.keys():
    Prob[(i,j)]=(a-b)*(Prob[(i,j)]-B)/(A-B) + b
print('Probability dictionary created')
#Create R
D={}
for (i,j) in Prob.keys():
    p=Prob[(i,j)]
    q=0.5*p
    D[(i,j)]=(-np.log(p),-np.log(q)+np.log(p))


print('D created')      
vdim=len(D)
R={}

for (i,j) in D.keys():
    (l,r)=D[(i,j)]
    k=0
    for (u,v) in D.keys():
        if i==u and j==v:
            R[(i,j,k)]=r
        else:
            R[(i,j,k)]=0
        k=k+1

print('R created')
B=25 # One percent of the nodes
Prob=None

opt=SolverFactory('gurobi')
opt.options['TimeLimit']=3600

#%%
times=[]

#randomly sample sources and sinks

for loop in range(1,11):
    print(f'Loop {loop}:')
    S=random.choice(list(N))
    T=random.choice(list(N-{S}))
    
    
    start=time.time()
    M=create_asymmetric_uncertainty_shortest_path_interdiction(D,N,R,S,T,vdim,B)
    end=time.time()
    build_time=end-start
    print(f'Finished creating model: {build_time} seconds')
    start=time.time()
    opt.solve(M,tee=True)
    end=time.time()
    times.append(end-start)
    print(f'Finished solving model: {end-start} seconds')
    


#%%    
k=0
times_array=np.empty((10,1))
for i in times:
    times_array[k]=i    
    k+=1
fig1, ax1 = plt.subplots()
ax1.set_title('Solve times')
ax1.boxplot(times_array)    
    
    
    
    
    

