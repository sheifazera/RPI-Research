import numpy as np
import itertools as it
import matplotlib.pyplot as plt
import random
import time
from pyomo.environ import *
from prettytable import PrettyTable, ALL
from shortestpath import create_asymmetric_uncertainty_shortest_path_interdiction
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
    D[(i,j)]=(-np.log(p),-np.log(q)+np.log(p))
      
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


B=25 #approximately 10 percent of the arcs can be interdicted

I=set(it.combinations(N,2))
times=[]

table=PrettyTable()
table.hrules=ALL
table.field_names=["Size of I", "Solve Time"]
loop=0
opt=SolverFactory('gurobi')
opt.options['TimeLimit']=1200
for (S,T) in I:
    print(f'Starting loop {loop} for ({S},{T})')
    M=create_asymmetric_uncertainty_shortest_path_interdiction(D,N,R,S,T,vdim,B)
    start=time.time()
    opt.solve(M)
    end=time.time()
    times.append(end-start)
    table.add_row([(S,T),end-start])
    loop=loop+1
    
    
#%%

k=0
times_array=np.empty((2701,1))
for i in times:
    times_array[k]=i    
    k+=1
fig1, ax1 = plt.subplots()
ax1.set_title('Solve times')
ax1.boxplot(times_array,showfliers=False)

print(table)

table_txt=table.get_string()
with open('Eastern_mass_ST_tables.txt','a') as file:
    file.write(table_txt)
    
           
    
    
    



