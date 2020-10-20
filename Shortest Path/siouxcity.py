#Sioux Falls Network
import numpy as np
from pyomo.environ import *
from shortestpath import create_shortest_path, create_shortest_path_robust_ellipsoid, return_path #, return_path_robust
matrix = np.loadtxt('siouxfalls.txt', usecols=(0,1,3))
#print(matrix) 

D = {(int(a),int(b)) : c for a,b,c in matrix}
N=set(range(1,25))
s=3
t=20
#print(D)

# %%
M=create_shortest_path(D,N,s,t)
opt=SolverFactory('gurobi')
results=opt.solve(M)
M.pprint()
path=return_path(M,D,N,s,t)
#M.pprint()

# %% Robust Sioux Falls
udim=7
N_disrupt={8,10,11,15,16,20,22}
P={}
i=0
for n in N_disrupt: #for all nodes disrupted (uncertainty dimension)
    for (u,v) in D.keys(): #check all edges 
        if u==n or v==n: #edge enters the disrupted node 
            P[(u,v,i)]=0.5*D[(u,v)] 
        elif (u,n) in D.keys() or (v,n) in D.keys(): #u or v is adjacent to the disrupted node
            P[(u,v,i)]=0.25*D[(u,v)]
        else:
            P[(u,v,i)]=0
    i=i+1
ME=create_shortest_path_robust_ellipsoid(D,N,P,s,t,udim)
results_robust_ellipsoid=opt.solve(ME) #,tee=True)
path_robust=return_path_robust(ME,D,N,P,s,t,udim)
ME.pprint()


#%% Robust Sioux Falls with Certain Interdiction Addition
B=2

udim=7
D = {(int(a),int(b)) : (c, 1000) for a,b,c in matrix}
N_disrupt={8,10,11,15,16,20,22}
P={}
i=0
for n in N_disrupt: #for all nodes disrupted (uncertainty dimension)
    for (u,v) in D.keys(): #check all edges 
        (l,r)=D[(u,v)]
        if u==n or v==n: #edge enters the disrupted node 
            P[(u,v,i)]=0.5*l 
        elif (u,n) in D.keys() or (v,n) in D.keys(): #u or v is adjacent to the disrupted node
            P[(u,v,i)]=0.25*l
        else:
            P[(u,v,i)]=0
    i=i+1
MI=create_shortest_path_robust_ellipsoid_certain_interdiction(D,N,P,s,t,udim, B)
opt=SolverFactory('gurobi')
results_interdiction=opt.solve(MI) #,tee=True)
path_interdiction=return_path_interdiction_robust(MI,D,N,P,s,t,udim,opt)
#MI.pprint()
print('Robust Sioux Falls with Certain Interdiction Additions')
MI.x.pprint()
print(path_interdiction)

# %% Robust Sioux Falls with Uncertain Interdiction Additions
B=2

udim=7
D = {(int(a),int(b)) : (c, 0.5*c) for a,b,c in matrix}
N_disrupt={8,10,11,15,16,20,22}
P={}
i=0
for n in N_disrupt: #for all nodes disrupted (uncertainty dimension)
    for (u,v) in D.keys(): #check all edges 
        (l,r)=D[(u,v)]
        if u==n or v==n: #edge enters the disrupted node 
            P[(u,v,i)]=0.5*l
        elif (u,n) in D.keys() or (v,n) in D.keys(): #u or v is adjacent to the disrupted node
            P[(u,v,i)]=0.25*l
        else:
            P[(u,v,i)]=0
    i=i+1
vdim=3    
R={}
for i in range(0,vdim):
    for (u,v) in D.keys(): #assign a value for all edges
        if i==0:
            if u in {1,2,3,4,5,6}:
                R[(u,v,i)]=2
            else:
                R[(u,v,i)]=0
        elif i==1:
            if u in {7,8,9,20,11,12,16, 17, 18}:
                R[(u,v,i)]=2
            else:
                R[(u,v,i)]=0
        elif i==2:
            if u in {13,14,15,19,20,21,22,23,24}:
                R[(u,v,i)]=2
            else:
                R[(u,v,i)]=0

    
    
MIU=create_shortest_path_robust_ellipsoid_uncertain_interdiction(D,N,P,R,s,t,udim,vdim, B)
opt=SolverFactory('gurobi')
results_interdiction_uncertain=opt.solve(MIU) #,tee=True)
path_interdiction_uncertain=return_path_interdiction_robust_uncertain(MIU,D,N,P,R,s,t,udim,vdim,opt)
#MIU.pprint()
print('Robust Sioux Falls with Certain Interdiction Additions')
MIU.x.pprint()
MIU.v.pprint()
print(path_interdiction_uncertain)




