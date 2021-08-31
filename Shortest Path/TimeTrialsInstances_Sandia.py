# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 10:17:16 2021

@author: sheifazera
"""

import numpy as np
import itertools as it
import random
import time
from pyomo.environ import *
#from prettytable import PrettyTable, ALL
from shortestpath_networkx import create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST
#%%
matrix = np.loadtxt('eastern_mass.txt', usecols=(0,1,3))
#print(matrix) 

Prob = {(int(a),int(b)) : c for a,b,c in matrix}
N=set(range(1,75)) #check this

A=max(Prob.values())
B=min(Prob.values())
a=1
b=0.1

for (i,j) in Prob.keys():
    Prob[(i,j)]=(a-b)*(Prob[(i,j)]-B)/(A-B) + b


#Create R
D={}
for (i,j) in Prob.keys():
    p=Prob[(i,j)]
    q=0.5*p
    G.add_edge(i,j,length=-np.log(p),interdicted_length=-np.log(q)+np.log(p))
      
vdim=len(D)
R={}

for (i,j) in G.edges:
    (l,r)=D[(i,j)]
    k=0
    for (u,v) in D.keys():
        if i==u and j==v:
            R[(i,j,k)]=r
        else:
            R[(i,j,k)]=0
        k=k+1


B=25 #approximately 10 percent of the arcs can be interdicted


random.seed(a=631996)
    
#Sample Instances
T=22 #Boston
time_1=set()
time_2=set()
time_4=set()
time_8=set()
time_16=set()
time_32=set()
#time_64=set()


opt=SolverFactory('gurobi_direct')
opt.options['TimeLimit']=1200

#%%
for loop in range(1,21):
#Add solves and timing somewhere
    #table=PrettyTable()
    #table.hrules=ALL
    #table.field_names=["Size of I", "Model Build Time", "Solve Time"]
    print(f'Starting loop {loop}')

    S_64=random.sample(list(N-{T}),64)
    '''
    I_64=set()
    for i in S_64:
        I_64.add((i,T))
    
    start=time.time()
    M=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST(D,N,R,I_64,vdim,B)
    end=time.time()
    build_time=end-start
    print(f'Finished creating model for 64 (S,T) pairs: {build_time} seconds')
    start=time.time()
    opt.solve(M,tee=True)
    end=time.time()
    time_64.add(end-start)
    print(f'Finished solving model for 64 (S,T) pairs: {end-start} seconds')
    table.add_row([64,build_time,end-start])
    '''
    S_32=random.sample(list(S_64), 32)
    I_32=set()
    for i in S_32:
        I_32.add((i,T))
    
    start=time.time()
    M=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST(D,N,R,I_32,vdim,B)
    end=time.time()
    build_time=end-start
    print(f'Finished creating model for 32 (S,T) pairs: {build_time} seconds')
    start=time.time()
    opt.solve(M,tee=True)
    end=time.time()
    time_32.add(end-start)
    print(f'Finished solving model for 32 (S,T) pairs: {end-start} seconds')
    #table.add_row([32,build_time,end-start])
    
    
    S_16=random.sample(list(S_32),16)
    I_16=set()
    for i in S_16:
        I_16.add((i,T))
    
    start=time.time()
    M=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST(D,N,R,I_16,vdim,B)
    end=time.time()
    build_time=end-start
    print(f'Finished creating model for 16 (S,T) pairs: {build_time} seconds')
    start=time.time()
    opt.solve(M)
    end=time.time()
    time_16.add(end-start)
    print(f'Finished solving model for 16 (S,T) pairs: {end-start} seconds')
    #table.add_row([16,build_time,end-start])
    
    
    S_8=random.sample(list(S_16),8)
    I_8=set()
    for i in S_8:
        I_8.add((i,T))
    
    start=time.time()
    M=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST(D,N,R,I_8,vdim,B)
    end=time.time()
    build_time=end-start
    print(f'Finished creating model for 8 (S,T) pairs: {build_time} seconds')
    start=time.time()
    opt.solve(M)
    end=time.time()
    time_8.add(end-start)
    print(f'Finished solving model for 8 (S,T) pairs: {end-start} seconds')
    #table.add_row([8,build_time,end-start])
    
    S_4=random.sample(list(S_8),4)
    I_4=set()
    for i in S_4:
        I_4.add((i,T))
    
    
    start=time.time()
    M=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST(D,N,R,I_4,vdim,B)
    end=time.time()
    build_time=end-start
    print(f'Finished creating model for 4 (S,T) pairs: {build_time} seconds')
    start=time.time()
    opt.solve(M)
    end=time.time()
    time_4.add(end-start)
    print(f'Finished solving model for 4 (S,T) pairs: {end-start} seconds')
    #table.add_row([4,build_time,end-start])
    
    S_2=random.sample(list(S_4),2)
    I_2=set()
    for i in S_2:
        I_2.add((i,T))
    
    start=time.time()
    M=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST(D,N,R,I_2,vdim,B)
    end=time.time()
    build_time=end-start
    print(f'Finished creating model for 2 (S,T) pairs: {build_time} seconds')
    start=time.time()
    opt.solve(M)
    end=time.time()
    time_2.add(end-start)
    print(f'Finished solving model for 2 (S,T) pairs: {end-start} seconds')
    #table.add_row([2,build_time,end-start])
    
    
    S_1=random.sample(list(S_2),1)       
    I_1=set()
    for i in S_1:           
        I_1.add((i,T))             
    
    start=time.time()
    M=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST(D,N,R,I_1,vdim,B)
    end=time.time()
    build_time=end-start
    print(f'Finished creating model for 1 (S,T) pair: {build_time} seconds')
    start=time.time()
    opt.solve(M)
    end=time.time()
    time_1.add(end-start)
    print(f'Finished solving model for 1 (S,T) pair: {end-start} seconds')
    #table.add_row([1,build_time,end-start])
    
    
    #table_txt=table.get_string()
    #with open('TimeTrials_MST_tables.txt','a') as file:
        #file.write('New Loop \n')
        #file.write(table_txt)
    
    #print(table)
    #Set up some sort of pretty table results
    
#%%
times=np.empty((20,5))
k=0 
for i in time_1:
    times[k,0]=i
    k+=1
k=0
for i in time_2:
    times[k,1]=i
    k+=1
k=0
for i in time_4:
    times[k,2]=i    
    k+=1
k=0
for i in time_8:
    times[k,3]=i  
    k+=1
k=0
for i in time_16:
    times[k,4]=i   
    k+=1
k=0

times_32_array=np.empty((20,1))

for i in time_32:
    times_32_array[k]=i    
    k+=1
    
set_label=np.array([1,2,4,8,16])
label_32=np.array([32])

import matplotlib.pyplot as plt  

fig1, ax1 = plt.subplots()
ax1.set_title('Solve times by number of (S,T) instances')
ax1.boxplot(times,labels=set_label)

fig2, ax2 = plt.subplots()
ax2.set_title('Solve times for 32 (S,T) pairs')
ax2.boxplot(times_32_array,labels=label_32)



#%% Not Going to Boston

matrix = np.loadtxt('eastern_mass.txt', usecols=(0,1,3))
#print(matrix) 

Prob = {(int(a),int(b)) : c for a,b,c in matrix}
N=set(range(1,75)) #check this

A=max(Prob.values())
B=min(Prob.values())
a=1
b=0.1

for (i,j) in Prob.keys():
    Prob[(i,j)]=(a-b)*(Prob[(i,j)]-B)/(A-B) + b


#Create R
for (i,j) in Prob.keys():
    p=Prob[(i,j)]
    q=0.5*p
    G.add_edge(i,j,length=-np.log(p),interdicted_length=-np.log(q)+np.log(p))
      
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

B=5 #approximately 10 percent of the arcs can be interdicted


random.seed(a=631996)
    
#Sample Instances
time_1=set()
time_2=set()
time_4=set()
time_8=set()
time_16=set()
time_32=set()
#time_64=set()


opt=SolverFactory('gurobi_direct')
opt.options['TimeLimit']=1200

ST=set(it.permutations(N,2))

#%%
for loop in range(1,21):
#Add solves and timing somewhere
    #table=PrettyTable()
    #table.hrules=ALL
    #table.field_names=["Size of I", "Model Build Time", "Solve Time"]
    print(f'Starting loop {loop}')

    #S_64=random.sample(list(N-{T}),64)
    
    I_64=set(random.sample(ST,64))
    '''
    for (S,T) in I_64:
        print(f'({S},{T})')
    '''
    '''
    I_64=set()
    for i in S_64:
        I_64.add((i,T))
    
    start=time.time()
    M=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST(D,N,R,I_64,vdim,B)
    end=time.time()
    build_time=end-start
    print(f'Finished creating model for 64 (S,T) pairs: {build_time} seconds')
    start=time.time()
    opt.solve(M)
    end=time.time()
    time_64.add(end-start)
    print(f'Finished solving model for 64 (S,T) pairs: {end-start} seconds')
    table.add_row([64,build_time,end-start])
    '''
    
    I_32=set(random.sample(list(I_64), 32))
    '''
    for (S,T) in I_32:
        print(f'({S},{T})')
    
    I_32=set()
    for i in S_32:
        I_32.add((i,T))
    '''
    start=time.time()
    M=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_32,vdim,B)
    end=time.time()
    build_time=end-start
    print(f'Finished creating model for 32 (S,T) pairs: {build_time} seconds')
    start=time.time()
    opt.solve(M)
    end=time.time()
    time_32.add(end-start)
    print(f'Finished solving model for 32 (S,T) pairs: {end-start} seconds')
    #table.add_row([32,build_time,end-start])
    
    
    I_16=set(random.sample(list(I_32),16))
    '''
    for (S,T) in I_16:
        print(f'({S},{T})')
    
    I_16=set()
    for i in S_16:
        I_16.add((i,T))
    '''
    start=time.time()
    M=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_16,vdim,B)
    end=time.time()
    build_time=end-start
    print(f'Finished creating model for 16 (S,T) pairs: {build_time} seconds')
    start=time.time()
    opt.solve(M)
    end=time.time()
    time_16.add(end-start)
    print(f'Finished solving model for 16 (S,T) pairs: {end-start} seconds')
    #table.add_row([16,build_time,end-start])
    
    
    I_8=set(random.sample(list(I_16),8))
    '''
    for (S,T) in I_8:
        print(f'({S},{T})')
    
    I_8=set()
    for i in S_8:
        I_8.add((i,T))
    '''
    start=time.time()
    M=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_8,vdim,B)
    end=time.time()
    build_time=end-start
    print(f'Finished creating model for 8 (S,T) pairs: {build_time} seconds')
    start=time.time()
    opt.solve(M)
    end=time.time()
    time_8.add(end-start)
    print(f'Finished solving model for 8 (S,T) pairs: {end-start} seconds')
    #table.add_row([8,build_time,end-start])
    
    I_4=set(random.sample(list(I_8),4))
    '''
    for (S,T) in I_4:
        print(f'({S},{T})')
    
    I_4=set()
    for i in S_4:
        I_4.add((i,T))
    '''
    
    start=time.time()
    M=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_4,vdim,B)
    end=time.time()
    build_time=end-start
    print(f'Finished creating model for 4 (S,T) pairs: {build_time} seconds')
    start=time.time()
    opt.solve(M)
    end=time.time()
    time_4.add(end-start)
    print(f'Finished solving model for 4 (S,T) pairs: {end-start} seconds')
    #table.add_row([4,build_time,end-start])
    
    I_2=set(random.sample(list(I_4),2))
    '''
    for (S,T) in I_2:
        print(f'({S},{T})')
    
    I_2=set()
    for i in S_2:
        I_2.add((i,T))
    '''
    start=time.time()
    M=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_2,vdim,B)
    end=time.time()
    build_time=end-start
    print(f'Finished creating model for 2 (S,T) pairs: {build_time} seconds')
    start=time.time()
    opt.solve(M)
    end=time.time()
    time_2.add(end-start)
    print(f'Finished solving model for 2 (S,T) pairs: {end-start} seconds')
    #table.add_row([2,build_time,end-start])
    
    
    I_1=set(random.sample(list(I_2),1))
    '''
    for (S,T) in I_1:
        print(f'({S},{T})')
         
    I_1=set()
    for i in S_1:           
        I_1.add((i,T))             
        '''    
    start=time.time()
    M=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I_1,vdim,B)
    end=time.time()
    build_time=end-start
    print(f'Finished creating model for 1 (S,T) pair: {build_time} seconds')
    start=time.time()
    opt.solve(M)
    end=time.time()
    time_1.add(end-start)
    print(f'Finished solving model for 1 (S,T) pair: {end-start} seconds')
    #table.add_row([1,build_time,end-start])
    
    
    #table_txt=table.get_string()
    #with open('TimeTrials_MST_tables.txt','a') as file:
    #   file.write('New Loop \n')
    #    file.write(table_txt)
    
    #print(table)
    #Set up some sort of pretty table results
    
#%%
times=np.empty((27,5))
k=0 
for i in time_1:
    times[k,0]=i
    k+=1
k=0
for i in time_2:
    times[k,1]=i
    k+=1
k=0
for i in time_4:
    times[k,2]=i    
    k+=1
k=0
for i in time_8:
    times[k,3]=i  
    k+=1
k=0
for i in time_16:
    times[k,4]=i   
    k+=1
k=0

times_32_array=np.empty((27,1))

for i in time_32:
    times_32_array[k]=i    
    k+=1
    
set_label=np.array([1,2,4,8,16])
label_32=np.array([32])

import matplotlib.pyplot as plt  

fig1, ax1 = plt.subplots()
ax1.set_title('Solve times by number of (S,T) instances, B=5')
ax1.boxplot(times,labels=set_label,showfliers=False)

fig2, ax2 = plt.subplots()
ax2.set_title('Solve times for 32 (S,T) pairs, B=5')
ax2.boxplot(times_32_array,labels=label_32, showfliers=False)


#%%

opt=SolverFactory('gurobi')
opt.options['TimeLimit']=1200
times=[]
ST=set(it.permutations(N,2))

for loop in range(1,21):
    I=set(random.sample(list(ST),24))
    start=time.time()
    M=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST(D,N,R,I,vdim,B)
    end=time.time()
    build_time=end-start
    print(f'Finished creating model for 24 (S,T) pairs: {build_time} seconds')
    start=time.time()
    opt.solve(M)
    end=time.time()
    times.append(end-start)
    print(f'Finished solving model for 24 (S,T) pairs: {end-start} seconds')

time_array=np.empty((20,1))
    #%%
k=0
for i in times:
    time_array[k]=i    
    k+=1
    
fig1, ax1 = plt.subplots()
ax1.set_title('Solve times for 24 (S,T) pairs')
ax1.boxplot(time_array)