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



# %% Asymmetric Uncertainty  


D = {(int(a),int(b)) : (c, c) for a,b,c in matrix}
B = 20
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

s=3
t=20

M=create_asymmetric_uncertainty_shortest_path_interdiction(D,N,R,s,t,vdim,B)
opt=SolverFactory('gurobi')
results=opt.solve(M)
print('Asymmetric Uncertainty Interdiction Results')
print('Interdiction Variable')
M.x.pprint()
print('Evader Path Variable')
M.w.pprint()
M.z.pprint()
print('Uncertainty Variable')
M.t.pprint()
path=return_path_asymmetric(M,D,N,s,t)
print(path)
print(f'Objective={value(M.Obj)}')
print(f'd_t={M.d[t].value}')


# %% Probability-based Sioux Falls
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
M=create_shortest_path_interdiction(D,N,s,t, B)   
opt=SolverFactory('gurobi')
results=opt.solve(M) 
path=return_path_interdiction(M,D,N,s,t,opt)    
print('Interdiction')
print('Interdiction Variable')
M.x.pprint()
M.d.pprint()
print(path)
print(f'Probability of Evasion={exp(-M.d[t].value)}')    


# %% Robust Super Source Sioux Falls

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
M=create_asymmetric_uncertainty_shortest_path_interdiction(D,N,R,s,t,vdim,B)
opt=SolverFactory('gurobi')
results=opt.solve(M)
print('Asymmetric Uncertainty Interdiction Results')
print('Interdiction Variable')
M.x.pprint()
print('Evader Path Variable')
M.w.pprint()
M.z.pprint()
#print('Uncertainty Variable')
#M.t.pprint()
path=return_path_asymmetric(M,D,N,s,t)
print(path)
print(f'Objective={value(M.Obj)}')
print(f'Nominal Probability of Evasion={exp(-M.d[t].value)}')      

# %% Multiple (S,T) Pairs
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

M=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST(D,N,R,I,vdim,B)
opt=SolverFactory('gurobi')
results=opt.solve(M) 
#path=return_path_interdiction(M,D,N,s,t,opt)    
print('Interdiction')
print('Interdiction Variable')
M.x.pprint()
#print(path)
for (S,T) in I:
    path=return_path_asymmetric_multiple_ST(M,D,N,S,T)
    print(f'Probability of Evasion on ({S},{T})={exp(-M.STBlocks[(S,T)].d[T].value)}')    
    print(path)



# %% Robust Probability-based Sioux Falls 
    

        