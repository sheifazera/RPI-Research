# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 10:49:44 2021

@author: sheifazera
"""

import numpy as np
import time
from prettytable import PrettyTable, ALL
import networkx as nx
from pyomo.environ import *
from shortestpath import create_shortest_path, create_shortest_path_robust_ellipsoid, return_path , create_asymmetric_uncertainty_shortest_path_interdiction #, return_path_robust


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
S=11
T=560

opt = SolverFactory('gurobi_persistent')
Budgets={10, 50, 100, 200, 500, 1000, 1500}

table=PrettyTable()
table.hrules=ALL
table.field_names=["Number of Arcs", "Budget", "Time to Create Model", "Gurobi Solve Time", "Total Solve Time"]

#%% 500 Arcs

# Randomly select 500 arcs for 132-ish nodes 
number_of_rows = Prob_array.shape[0]
random_indices = np.random.choice(number_of_rows, size=500, replace=False)
Prob_500_array=Prob_array[random_indices, :]

Prob_500={}
for k in range(0,len(Prob_500_array)):
    (i,j)=Prob_500_array[k,0]
    Prob_500[(i,j)]=Prob_500_array[k,1]
    
    
# Use Networkx to check that there IS a path from source to sink
DG = nx.DiGraph()
for (i,j) in Prob_500.keys():
    DG.add_weighted_edges_from([(i, j, Prob_500[(i,j)])])
flag=nx.has_path(DG, S, T) 
if flag==False:
    raise Exception('There is no path from the source to the sink for the selected subset of arcs')
 
# Use my method to solve this
D_500={}
for (i,j) in Prob_500.keys():
    p=Prob_500[(i,j)]
    q=0.5*p
    D_500[(i,j)]=(-np.log(p),-np.log(q)+np.log(p))
    

vdim_500=len(D_500)
R_500={}

for (i,j) in D_500.keys():
    (l,r)=D_500[(i,j)]
    k=0
    for (u,v) in D_500.keys():
        if i==u and j==v:
            R_500[(i,j,k)]=r
        else:
            R_500[(i,j,k)]=0
        k=k+1
             
print('Complete: 500 Arc Dictionaries')
Budgets_500={10,50, 100,200}
for B in Budgets_500:
    start = time.time()
    M_500=create_asymmetric_uncertainty_shortest_path_interdiction(D_500,N,R_500,S,T,vdim_500,B)
    end = time.time()
    model_time_500=end-start
    opt.options['TimeLimit']=600
    opt.set_instance(M_500)
    start = time.time()
    results_500=opt.solve()
    end = time.time()
    total_time_500=end-start
    solve_time_500=results_500.solver.wallclock_time  
    table.add_row([500, B, model_time_500, solve_time_500,total_time_500])
    print(f'Finished Solving at budget={B}')
    
M_500=None
D_500=None
R_500=None

#%% 1000 Arcs 

# Randomly select 1000 arcs for 132-ish nodes 
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
flag=nx.has_path(DG, S, T) 
if flag==False:
    raise Exception('There is no path from the source to the sink for the selected subset of arcs')
 
# Use my method to solve this
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
'''             
B=10 

M_1000=create_asymmetric_uncertainty_shortest_path_interdiction(D,N,R,S,T,vdim,B)
start = time.time()
results_1000=opt.solve(M_1000)
end = time.time()
time_1000=end-start
print(end - start)
table.add_row([1000,time_1000]) 
 '''        
print('Complete: 1000 Arc Dictionaries')
Budgets_1000={10,50, 100,200}
for B in Budgets_1000:
    start = time.time()
    M_1000=create_asymmetric_uncertainty_shortest_path_interdiction(D_1000,N,R_1000,S,T,vdim_1000,B)
    end = time.time()
    model_time_1000=end-start
    opt = SolverFactory('gurobi_persistent')
    opt.options['TimeLimit']=600
    opt.set_instance(M_1000)
    start = time.time()
    results_1000=opt.solve()
    end = time.time()
    total_time_1000=end-start
    solve_time_1000=results_1000.solver.wallclock_time  
    table.add_row([1000, B, model_time_1000, solve_time_1000,total_time_1000])
    print(f'Finished Solving at budget={B}')

M_1000=None
D_1000=None
R_1000=None
#%% 5000 Arcs


# Randomly select 5000 arcs for 132-ish nodes 
number_of_rows = Prob_array.shape[0]
random_indices = np.random.choice(number_of_rows, size=5000, replace=False)
Prob_5000_array=Prob_array[random_indices, :]

Prob_5000={}
for k in range(0,len(Prob_5000_array)):
    (i,j)=Prob_5000_array[k,0]
    Prob_5000[(i,j)]=Prob_5000_array[k,1]
    
    
# Use Networkx to check that there IS a path from source to sink
DG = nx.DiGraph()
for (i,j) in Prob_5000.keys():
    DG.add_weighted_edges_from([(i, j, Prob_5000[(i,j)])])
flag=nx.has_path(DG, S, T) 
if flag==False:
    raise Exception('There is no path from the source to the sink for the selected subset of arcs')
 
# Use my method to solve this
D_5000={}
for (i,j) in Prob_5000.keys():
    p=Prob_5000[(i,j)]
    q=0.5*p
    D_5000[(i,j)]=(-np.log(p),-np.log(q)+np.log(p))
    

vdim_5000=len(D_5000)
R_5000={}

for (i,j) in D_5000.keys():
    (l,r)=D_5000[(i,j)]
    k=0
    for (u,v) in D_5000.keys():
        if i==u and j==v:
            R_5000[(i,j,k)]=r
        else:
            R_5000[(i,j,k)]=0
        k=k+1
'''            
B=10 

M_5000=create_asymmetric_uncertainty_shortest_path_interdiction(D,N,R,S,T,vdim,B)
start = time.time()
results_5000=opt.solve(M_5000)
end = time.time()
time_5000=end-start
print(end - start)  
table.add_row([5000,time_5000])   
'''
print('Complete: 5000 Arc Dictionaries')
for B in Budgets:
    start = time.time()
    M_5000=create_asymmetric_uncertainty_shortest_path_interdiction(D_5000,N,R_5000,S,T,vdim_5000,B)
    end = time.time()
    model_time_5000=end-start
    opt = SolverFactory('gurobi_persistent')
    opt.options['TimeLimit']=600
    opt.set_instance(M_5000)
    start = time.time()
    results_5000=opt.solve()
    end = time.time()
    total_time_5000=end-start
    solve_time_5000=results_5000.solver.wallclock_time  
    table.add_row([500, B, model_time_5000, solve_time_5000,total_time_5000])
    print(f'Finished Solving at budget={B}')

M_5000=None
D_5000=None
R_5000=None
#%% 10000 Arcs

# Randomly select 10000 arcs for 132-ish nodes 
number_of_rows = Prob_array.shape[0]
random_indices = np.random.choice(number_of_rows, size=10000, replace=False)
Prob_10000_array=Prob_array[random_indices, :]

Prob_10000={}
for k in range(0,len(Prob_10000_array)):
    (i,j)=Prob_10000_array[k,0]
    Prob_10000[(i,j)]=Prob_10000_array[k,1]
    
    
# Use Networkx to check that there IS a path from source to sink
DG = nx.DiGraph()
for (i,j) in Prob_10000.keys():
    DG.add_weighted_edges_from([(i, j, Prob_10000[(i,j)])])
flag=nx.has_path(DG, S, T) 
if flag==False:
    raise Exception('There is no path from the source to the sink for the selected subset of arcs')
 
# Use my method to solve this
D_10000={}
for (i,j) in Prob_10000.keys():
    p=Prob_10000[(i,j)]
    q=0.5*p
    D_10000[(i,j)]=(-np.log(p),-np.log(q)+np.log(p))
    

vdim_10000=len(D_10000)
R_10000={}

for (i,j) in D_10000.keys():
    (l,r)=D_10000[(i,j)]
    k=0
    for (u,v) in D_10000.keys():
        if i==u and j==v:
            R_10000[(i,j,k)]=r
        else:
            R_10000[(i,j,k)]=0
        k=k+1
'''             
B=10 

M_10000=create_asymmetric_uncertainty_shortest_path_interdiction(D,N,R,S,T,vdim,B)
start = time.time()
results_10000=opt.solve(M_10000)
end = time.time()
time_10000=end-start
print(end - start) 
table.add_row([10000,time_10000]) 
'''
print('Complete: 10000 Arc Dictionaries')
for B in Budgets:
    start = time.time()
    M_10000=create_asymmetric_uncertainty_shortest_path_interdiction(D_10000,N,R_10000,S,T,vdim_10000,B)
    end = time.time()
    model_time_10000=end-start
    opt = SolverFactory('gurobi_persistent')
    opt.options['TimeLimit']=600
    opt.set_instance(M_10000)
    start = time.time()
    results_10000=opt.solve()
    end = time.time()
    total_time_10000=end-start
    solve_time_10000=results_10000.solver.wallclock_time  
    table.add_row([10000, B, model_time_10000, solve_time_10000,total_time_10000])
    print(f'Finished Solving at budget={B}')

M_10000=None
D_10000=None
R_10000=None
#%% All Arcs 

 
# Use my method to solve this
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
'''            
B=10 

M=create_asymmetric_uncertainty_shortest_path_interdiction(D,N,R,S,T,vdim,B)
start = time.time()
results=opt.solve(M)
end = time.time()
time=end-start
print(end - start) 
table.add_row([17338,time]) 
'''
print('Complete: Total Arc Dictionaries')
for B in Budgets:
    start = time.time()
    M=create_asymmetric_uncertainty_shortest_path_interdiction(D,N,R,S,T,vdim,B)
    end = time.time()
    model_time=end-start
    opt = SolverFactory('gurobi_persistent')
    opt.options['TimeLimit']=600
    opt.set_instance(M)
    start = time.time()
    results=opt.solve()
    end = time.time()
    total_time=end-start
    solve_time=results.solver.wallclock_time  
    table.add_row([17338, B, model_time, solve_time,total_time])
    print(f'Finished Solving at budget={B}')

M=None
D=None
R=None
#%%
print(table)
#%%
'''
print('Complete: Creating Dictionaries')

Budgets={10, 50, 100, 500, 1000, 1500}
for B in Budgets:
    table=PrettyTable()
    table.hrules=ALL
    table.field_names=["Number of Arcs", "Time to Create Model", "Gurobi Solve Time", "Total Solve Time"]
    if B <500:
        start = time.time()
        M_500=create_asymmetric_uncertainty_shortest_path_interdiction(D_500,N,R_500,S,T,vdim_500,B)
        end = time.time()
        model_time_500=end-start
        opt.set_instance(M_500)
        start = time.time()
        results_500=opt.solve()
        end = time.time()
        total_time_500=end-start
        solve_time_500=results_500.solver.wallclock_time  
        table.add_row([500,model_time_500, solve_time_500,total_time_500]) 
    else:
        table.add_row([500, "N/A", "N/A", "N/A"])
    
    print(f'Budget={B}, Number of Arcs=500')
    
    if B < 1500:
        start = time.time()
        M_1000=create_asymmetric_uncertainty_shortest_path_interdiction(D_1000,N,R_1000,S,T,vdim_1000,B)
        end = time.time()
        model_time_1000=end-start
        opt.set_instance(M_1000)
        start = time.time()
        results_1000=opt.solve()
        end = time.time()
        total_time_1000=end-start
        solve_time_1000=results_1000.solver.wallclock_time
        table.add_row([1000,model_time_1000, solve_time_1000,total_time_1000]) 
    else:
        table.add_row([1000,"N/A", "N/A", "N/A"])
        
    print(f'Budget={B}, Number of Arcs=1000')

    start = time.time()    
    M_5000=create_asymmetric_uncertainty_shortest_path_interdiction(D_5000,N,R_5000,S,T,vdim_5000,B)
    end = time.time()
    model_time_5000=end-start
    opt.set_instance(M_5000)
    start = time.time()
    results_5000=opt.solve()
    end = time.time()
    total_time_5000=end-start
    solve_time_5000=results_5000.solver.wallclock_time
    table.add_row([5000,model_time_5000, solve_time_5000,total_time_5000])

    print(f'Budget={B}, Number of Arcs=5000')
    start = time.time()
    M_10000=create_asymmetric_uncertainty_shortest_path_interdiction(D_10000,N,R_10000,S,T,vdim_10000,B)
    end = time.time()
    model_time_10000=end-start
    opt.set_instance(M_10000)
    start = time.time()
    results_10000=opt.solve()
    end = time.time()
    total_time_10000=end-start
    solve_time_10000=results_10000.solver.wallclock_time
    table.add_row([10000,model_time_10000, solve_time_10000,total_time_10000])

    print(f'Budget={B}, Number of Arcs=10000')
    start = time.time()
    M=create_asymmetric_uncertainty_shortest_path_interdiction(D,N,R,S,T,vdim,B)
    end = time.time()
    model_time=end-start
    opt.set_instance(M)
    start = time.time()
    results=opt.solve()
    end = time.time()
    total_time=end-start
    solve_time=results.solver.wallclock_time
    

    table.add_row([17338,model_time, solve_time,total_time]) 
    
    print(f'Budget={B}, Number of Arcs=Full')
    
    
    print(f'Budget={B}')
    print(table)
    '''