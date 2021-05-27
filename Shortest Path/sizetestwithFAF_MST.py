# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 10:29:42 2021

@author: sheifazera
"""

import numpy as np
import time
import random
from prettytable import PrettyTable, ALL
import networkx as nx
from pyomo.environ import *
from shortestpath import create_shortest_path, create_shortest_path_robust_ellipsoid, return_path , create_asymmetric_uncertainty_shortest_path_interdiction, create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST #, return_path_robust


#%% Input matrix 
matrix = np.genfromtxt('FAF_nodes.csv', delimiter=',', skip_header=1)

#%% DO NOT RELOAD THAT MONSTER AGAIN
np.random.seed(seed=631996) #To get the same rows everytime we run this
Prob = {(int(a),int(b)) : c for a,b,c in matrix}

N=set()
for (i,j) in Prob.keys():
    if i not in N:
        N.add(i)
    if j not in N:
        N.add(j)

for (i,j) in Prob.keys(): #Turn these weird values into probabilities
    if Prob[(i,j)] >1 and Prob[(i,j)] < 100:
        Prob[(i,j)]=Prob[(i,j)]/100 
    elif Prob[(i,j)] >= 100:
        Prob[(i,j)]=1
    elif Prob[(i,j)] <0.1:
        Prob[(i,j)]=0.1
        
data = list(Prob.items())
Prob_array = np.array(data, dtype=object)
    
matrix=None


#%%   
random.seed(a=631996)
T=560
S_2=set(random.choices(list(N),k=2))
S_4=set(random.choices(list(N),k=4))
S_8=set(random.choices(list(N),k=8))
S_16=set(random.choices(list(N),k=16))
S_16.remove(513)
S_16.add(511)
S_32=set(random.choices(list(N),k=32))
S_32.remove(261)
S_32.add(211)
S_64=set(random.choices(list(N),k=64))
S_64.remove(262)
S_64.add(230)


B=10
#%% 

number_of_rows = Prob_array.shape[0]
random_indices = np.random.choice(number_of_rows, size=1000, replace=False)
Prob_1000_array=Prob_array[random_indices, :]

Prob_1000={}
for k in range(0,len(Prob_1000_array)):
    (i,j)=Prob_1000_array[k,0]
    Prob_1000[(i,j)]=Prob_1000_array[k,1]
    
    
# Use Networkx to check that there IS a path from source to sink
DG = nx.DiGraph()
for (i,j) in Prob_1000.keys():
    DG.add_weighted_edges_from([(i, j, Prob_1000[(i,j)])])
#%% 
for S in S_2.union(S_4,S_8,S_16,S_32,S_64):
        flag=nx.has_path(DG, S, T) 
        if flag==False:
            raise Exception(f'There is no path from the source {S} to the sink for the selected subset of arcs')
#%%
I_2=set()
for i in S_2:
    I_2.add((i,T))
I_4=set()
for i in S_4:
    I_4.add((i,T))   
I_8=set()
for i in S_8:
    I_8.add((i,T))
I_16=set()
for i in S_16:
    I_16.add((i,T))
I_32=set()
for i in S_32:
    I_32.add((i,T))   
I_64=set()
for i in S_64:
    I_64.add((i,T))
#%% 

D_1000={}
for (i,j) in Prob_1000.keys():
    p=Prob_1000[(i,j)]
    q=0.5*p
    D_1000[(i,j)]=(-np.log(p),-np.log(q)+np.log(p))
    

vdim_1000=len(D_1000)
R_1000={}

for (i,j) in D_1000.keys():
    (l,r)=D_1000[(i,j)]
    k=0
    for (u,v) in D_1000.keys():
        if i==u and j==v:
            R_1000[(i,j,k)]=r
        else:
            R_1000[(i,j,k)]=0
        k=k+1

Prob_1000_array=None
Prob=None      
#%%

table=PrettyTable()
table.hrules=ALL
table.field_names=["Size of I",  "Gurobi Solve Time"]

opt=SolverFactory('gurobi_persistent')
opt.options['TimeLimit']=1200


M_2=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST(D_1000,N,R_1000,I_2,vdim_1000,B)
print('Finished Creating Model M_2')
opt.set_instance(M_2)
results=opt.solve()
time_2=results.solver.wallclock_time
table.add_row([2,time_2])
print('Finished Solving Model M_2')

M_4=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST(D_1000,N,R_1000,I_4,vdim_1000,B)
print('Finished Creating Model M_4')
opt.set_instance(M_4)
results=opt.solve()
time_4=results.solver.wallclock_time
table.add_row([4,time_4])
print('Finished Solving Model M_4')

M_8=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST(D_1000,N,R_1000,I_8,vdim_1000,B)
print('Finished Creating Model M_8')
opt.set_instance(M_8)
results=opt.solve()
time_8=results.solver.wallclock_time
table.add_row([8,time_8])
print('Finished Solving Model M_8')

M_16=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST(D_1000,N,R_1000,I_16,vdim_1000,B)  
print('Finished Creating Model M_16')  
opt.set_instance(M_16)
results=opt.solve()
time_16=results.solver.wallclock_time
table.add_row([16,time_16])
print('Finished Solving Model M_16')

M_32=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST(D_1000,N,R_1000,I_32,vdim_1000,B)
print('Finished Creating Model M_32')
opt.set_instance(M_32)
results=opt.solve()
time_32=results.solver.wallclock_time
table.add_row([32,time_32])
print('Finished Solving Model M_32')

M_64=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST(D_1000,N,R_1000,I_64,vdim_1000,B)
print('Finished Creating Model M_64')
opt.set_instance(M_64)
results=opt.solve()
time_64=results.solver.wallclock_time
table.add_row([64,time_64])
print('Finished Solving Model M_64')

#%%
print(table)

