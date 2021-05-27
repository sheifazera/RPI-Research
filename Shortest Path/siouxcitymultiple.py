# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 10:49:20 2021

@author: sheifazera


Sioux Falls with multiple sources and sinks, giving all paths in the result

"""

import numpy as np
import time
from pyomo.environ import *
from shortestpath import path_length, adjusted_dictionary, create_shortest_path, create_asymmetric_uncertainty_shortest_path_interdiction, create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST, create_shortest_path_robust_ellipsoid,return_path_interdiction, return_path, create_shortest_path_interdiction, find_paths, find_paths_MST #, return_path_robust
matrix = np.loadtxt('siouxfalls.txt', usecols=(0,1,3))
#print(matrix) 

D = {(int(a),int(b)) : c for a,b,c in matrix}
N=set(range(1,25))
s=3
t=20
#print(D)

from prettytable import PrettyTable, ALL

timetable=PrettyTable()
timetable.field_names=["Method", "Time to Run"]

# %% Shortest Path Interdiction with a Super Source
Prob = {(int(a),int(b)) : c for a,b,c in matrix}
for (i,j) in Prob.keys():
    if Prob[(i,j)] ==10:
        Prob[(i,j)] = (0.1,0.05) #p_ij, q_ij 
    elif Prob[(i,j)] ==9 :
        Prob[(i,j)] = (0.2, 0.1)
    elif Prob[(i,j)] ==8 :  
        Prob[(i,j)] = (0.3,0.15)
    elif Prob[(i,j)] ==7 :  
        Prob[(i,j)] = (0.4,0.2)
    elif Prob[(i,j)] ==6 :  
        Prob[(i,j)] = (0.5, 0.25)
    elif Prob[(i,j)] ==5 :  
        Prob[(i,j)] = (0.6, 0.3)
    elif Prob[(i,j)] ==4 :  
        Prob[(i,j)] = (0.7, 0.35)
    elif Prob[(i,j)] ==3 :  
        Prob[(i,j)] = (0.8, 0.4)
    elif Prob[(i,j)] ==2 :  
        Prob[(i,j)] = (0.9, 0.45)    
    else:
        raise Exception('Original Sioux Falls value out of range 2-10')



#Add the super source   100     
N.add(100)
N_external={1,2,3,7,12,13,18,20,21,24}
for i in N_external:
    Prob[(i,100)]=(1,1)
    Prob[(100,i)]=(1,1)


s=100
t=10
#Convert to lengths using \ell_ij=-\ln(p_ij)
D={}
for (i,j) in Prob.keys():
    (p,q)=Prob[(i,j)]
    D[(i,j)]=(-np.log(p),-np.log(q)+np.log(p))
    
B=5  
M_SS=create_shortest_path_interdiction(D,N,s,t, B)   

opt=SolverFactory('gurobi')
start=time.time()
results=opt.solve(M_SS)
end=time.time() 
timetable.add_row(['SS',end-start])
path_SS=return_path_interdiction(M_SS,D,N,s,t,opt)    
print('Interdiction')
print('Interdiction Variable')
tol = 1e-4
for k in M_SS.x.keys():
    if abs(M_SS.x[k].value)>tol:
        print(M_SS.x[k].getname(), M_SS.x[k].value)
M_SS.d.pprint()
print(path_SS)
print(f'Probability of Evasion={exp(-M_SS.d[t].value)}')    


# %% Robust Shortest Path Interdiction with a Super Source

#Uncertainty

vdim=len(D)
#vdim=0
R={}

for (i,j) in D.keys():
    (l,r)=D[(i,j)]
    k=0
    for (u,v) in D.keys():
        if i==u and j==v:
            if i==100 or j==100:
                R[(i,j,k)]=0 #The arcs from super source shouldn't be variable
            else:
                R[(i,j,k)]=r
        else:
            R[(i,j,k)]=0
        k=k+1


B=5
M_SSR=create_asymmetric_uncertainty_shortest_path_interdiction(D,N,R,s,t,vdim,B)
opt=SolverFactory('gurobi')
start=time.time()
results=opt.solve(M_SSR)
end=time.time() 
timetable.add_row(['SSR',end-start])
print('Asymmetric Uncertainty Interdiction Results')
print('Interdiction Variable')
tol = 1e-4
for k in M_SSR.x.keys():
    if abs(M_SSR.x[k].value)>tol:
        print(M_SSR.x[k].getname(), M_SSR.x[k].value)
print('Evader Path Variable')
for k in M_SSR.x.keys():
    if abs(M_SSR.w[k].value)>tol:
        print(M_SSR.w[k].getname(), M_SSR.w[k].value)
    if abs(M_SSR.z[k].value)>tol:
        print(M_SSR.z[k].getname(), M_SSR.z[k].value)
#print('Uncertainty Variable')
#M.t.pprint()
(no_paths_SSR,path_SSR)=find_paths(M_SSR,D,N,s,t)
for list in path_SSR:
    print(list)
print(f'Objective={value(M_SSR.Obj)}')
print(f'Nominal Probability of Evasion={exp(-M_SSR.d[t].value)}')
D_SSR=D      

# %% Shortest Path Interdiction with Multiple (S,T) Pairs
N=set(range(1,25))
Prob = {(int(a),int(b)) : c for a,b,c in matrix}
for (i,j) in Prob.keys():
    if Prob[(i,j)] ==10:
        Prob[(i,j)] = (0.1,0.05) #p_ij, q_ij 
    elif Prob[(i,j)] ==9 :
        Prob[(i,j)] = (0.2, 0.1)
    elif Prob[(i,j)] ==8 :  
        Prob[(i,j)] = (0.3,0.15)
    elif Prob[(i,j)] ==7 :  
        Prob[(i,j)] = (0.4,0.2)
    elif Prob[(i,j)] ==6 :  
        Prob[(i,j)] = (0.5, 0.25)
    elif Prob[(i,j)] ==5 :  
        Prob[(i,j)] = (0.6, 0.3)
    elif Prob[(i,j)] ==4 :  
        Prob[(i,j)] = (0.7, 0.35)
    elif Prob[(i,j)] ==3 :  
        Prob[(i,j)] = (0.8, 0.4)
    elif Prob[(i,j)] ==2 :  
        Prob[(i,j)] = (0.9, 0.45)    
    else:
        raise Exception('Original Sioux Falls value out of range 2-10')
   
N_external={1,2,3,7,12,13,18,20,21,24}
t=10
I=set()
for i in N_external:
    I.add((i,t))

#Convert to lengths using \ell_ij=-\ln(p_ij)
D={}
for (i,j) in Prob.keys():
    (p,q)=Prob[(i,j)]
    D[(i,j)]=(-np.log(p),-np.log(q)+np.log(p))
    
B=5  

vdim=0
R={}

M_MST=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST(D,N,R,I,vdim,B)
opt=SolverFactory('gurobi')
start=time.time()
results=opt.solve(M_MST) 
end=time.time() 
timetable.add_row(['MST',end-start])
#path=return_path_interdiction(M,D,N,s,t,opt)    
print('Interdiction')
print('Interdiction Variable')
tol = 1e-4
for k in M_MST.x.keys():
    if abs(M_MST.x[k].value)>tol:
        print(M_MST.x[k].getname(), M_MST.x[k].value)
#print(path)
for (S,T) in I:
    (no_paths_MST,path_MST)=find_paths_MST(M_MST,D,N,S,T)
    print(f'Probability of Evasion on ({S},{T})={exp(-M_MST.STBlocks[(S,T)].d[T].value)}')  
    for list in path_MST:
        print(list)
    
# %%  Robust Shortest Path Interdiction with Multiple (S,T) Pairs
Prob = {(int(a),int(b)) : c for a,b,c in matrix}
for (i,j) in Prob.keys():
    if Prob[(i,j)] ==10:
        Prob[(i,j)] = (0.1,0.05) #p_ij, q_ij 
    elif Prob[(i,j)] ==9 :
        Prob[(i,j)] = (0.2, 0.1)
    elif Prob[(i,j)] ==8 :  
        Prob[(i,j)] = (0.3,0.15)
    elif Prob[(i,j)] ==7 :  
        Prob[(i,j)] = (0.4,0.2)
    elif Prob[(i,j)] ==6 :  
        Prob[(i,j)] = (0.5, 0.25)
    elif Prob[(i,j)] ==5 :  
        Prob[(i,j)] = (0.6, 0.3)
    elif Prob[(i,j)] ==4 :  
        Prob[(i,j)] = (0.7, 0.35)
    elif Prob[(i,j)] ==3 :  
        Prob[(i,j)] = (0.8, 0.4)
    elif Prob[(i,j)] ==2 :  
        Prob[(i,j)] = (0.9, 0.45)    
    else:
        raise Exception('Original Sioux Falls value out of range 2-10')
   
N_external={1,2,3,7,12,13,18,20,21,24}
t=10
I=set()
for i in N_external:
    I.add((i,t))

#Convert to lengths using \ell_ij=-\ln(p_ij)
D={}
for (i,j) in Prob.keys():
    (p,q)=Prob[(i,j)]
    D[(i,j)]=(-np.log(p),-np.log(q)+np.log(p))
    
B=5  
vdim=len(D)
#vdim=0
R={}

for (i,j) in D.keys():
    (l,r)=D[(i,j)]
    k=0
    for (u,v) in D.keys():
        if i==u and j==v:
            if i==100 or j==100:
                R[(i,j,k)]=0 #The arcs from super source shouldn't be variable
            else:
                R[(i,j,k)]=r
        else:
            R[(i,j,k)]=0
        k=k+1

M_MSTR=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST(D,N,R,I,vdim,B)
opt=SolverFactory('gurobi')
start=time.time()
results=opt.solve(M_MSTR) 
end=time.time() 
timetable.add_row(['MSTR',end-start])
#path=return_path_interdiction(M,D,N,s,t,opt)    
print('Interdiction')
print('Interdiction Variable')
tol = 1e-4
for k in M_MSTR.x.keys():
    if abs(M_MSTR.x[k].value)>tol:
        print(M_MSTR.x[k].getname(), M_MSTR.x[k].value)
for (S,T) in I:
    (no_paths_MSTR, path)=find_paths_MST(M_MSTR,D,N,S,T)
    print(f'Probability of Evasion on ({S},{T})={exp(-M_MSTR.STBlocks[(S,T)].d[T].value)}')    
    for list in path:
        print(list)



# %% Getting a pretty table for comparison

from prettytable import PrettyTable, ALL

table=PrettyTable()
table.hrules=ALL
table.field_names=["(S,T)", "SS SP", "SS Prob", "SSR SP", "SSR Nom Prob", "SSR Robust Prob", "MST SP", "MST Prob", "MSTR SP", "MSTR Nom Prob", "MSTR Robust Prob"]

for (S,T) in I:
    
    #Shortest Path Interdiction with a Super Source (currently does not support multiple paths)
    if S in path_SS:
        SS1=path_SS
        SS2=exp(-M_SS.d[t].value)
    else:
         SS1="NA"
         SS2= "NA"
    
    # Robust Shortest Path Interdiction with a Super Source
    SSR1=[]
    SSR2=[]
    SSR3=[]
    for path in path_SSR:     
        if S in path: 
            path_str=[]
            for i in path:
                path_str.append(str(i))
            path_str="-".join(path_str)
            SSR1.append(path_str)
            SSR2.append(str(exp(-M_SSR.d[t].value)))
            D_adj=adjusted_dictionary(M_SSR,D_SSR,R,S,T)
            SSR3.append(str(exp(-path_length(M_SSR,path,D_adj,S,T))))
               
    SSR1="\n".join(SSR1)
    SSR2="\n".join(SSR2)
    SSR3="\n".join(SSR3) 
    (no_of_path_MST, path_MST)=find_paths_MST(M_MST,D,N,S,T)
    
    MST1=[]
    MST2=[]
    for path in path_MST:
        path_str=[]
        for i in path:
            path_str.append(str(i))
        path_str="-".join(path_str)
        MST1.append(path_str)
        
        MST2.append(str(exp(-M_MST.STBlocks[(S,T)].d[T].value)))
        
    MST1="\n".join(MST1) 
    MST2="\n".join(MST2)
    (no_of_path_MSTR, path_MSTR)=find_paths_MST(M_MSTR,D,N,S,T)   
    D_adj=adjusted_dictionary(M_MSTR,D,R,S,T)
    MSTR1=[]
    MSTR2=[]
    MSTR3=[]
    for path in path_MSTR:
        path_str=[]
        for i in path:
            path_str.append(str(i))
        path_str="-".join(path_str)
        MSTR1.append(path_str)
        
        MSTR2.append(str(exp(-M_MSTR.STBlocks[(S,T)].d[T].value)))
        
        
        evad=path_length(M_MSTR,path,D_adj,S,T)
        MSTR3.append(str(exp(-evad)))
    MSTR1="\n".join(MSTR1)    
    MSTR2="\n".join(MSTR2)    
    MSTR3="\n".join(MSTR3)    
    
    
    table.add_row([(S,T), SS1, SS2, SSR1, SSR2, SSR3,MST1,MST2, MSTR1,MSTR2, MSTR3])


print(f'Budget={B}')
print(table)

print(timetable)

