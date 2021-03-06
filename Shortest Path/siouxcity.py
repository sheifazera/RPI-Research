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
M_SS=create_shortest_path_interdiction(D,N,s,t, B)   
opt=SolverFactory('gurobi')
results=opt.solve(M_SS) 
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
M_SSR=create_asymmetric_uncertainty_shortest_path_interdiction(D,N,R,s,t,vdim,B)
opt=SolverFactory('gurobi')
results=opt.solve(M_SSR)
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
path_SSR=return_path_asymmetric(M_SSR,D,N,s,t)
print(path_SSR)
print(f'Objective={value(M_SSR.Obj)}')
print(f'Nominal Probability of Evasion={exp(-M_SSR.d[t].value)}')
D_SSR=D      

# %% Multiple (S,T) Pairs
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
results=opt.solve(M_MST) 
#path=return_path_interdiction(M,D,N,s,t,opt)    
print('Interdiction')
print('Interdiction Variable')
tol = 1e-4
for k in M_MST.x.keys():
    if abs(M_MST.x[k].value)>tol:
        print(M_MST.x[k].getname(), M_MST.x[k].value)
#print(path)
for (S,T) in I:
    path_MST=return_path_asymmetric_multiple_ST(M_MST,D,N,S,T)
    print(f'Probability of Evasion on ({S},{T})={exp(-M_MST.STBlocks[(S,T)].d[T].value)}')    
    print(path_MST)
    
# %%  Robust Multiple ST Pairs
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
results=opt.solve(M_MSTR) 
#path=return_path_interdiction(M,D,N,s,t,opt)    
print('Interdiction')
print('Interdiction Variable')
tol = 1e-4
for k in M_MSTR.x.keys():
    if abs(M_MSTR.x[k].value)>tol:
        print(M_MSTR.x[k].getname(), M_MSTR.x[k].value)

for (S,T) in I:
    path=return_path_asymmetric_multiple_ST(M_MSTR,D,N,S,T)
    print(f'Probability of Evasion on ({S},{T})={exp(-M_MSTR.STBlocks[(S,T)].d[T].value)}')    
    print(path)



# %% Getting a pretty table for comparison

from prettytable import PrettyTable

table=PrettyTable()
table.field_names=["(S,T)", "SS SP", "SS Prob", "SSR SP", "SSR Nom Prob", "SSR Robust Prob", "MST SP", "MST Prob", "MSTR SP", "MSTR Nom Prob", "MSTR Robust Prob"]
for (S,T) in I:
    if S in path_SS:
        SS1=path_SS
        SS2=exp(-M_SS.d[t].value)
    else:
         SS1="NA"
         SS2= "NA"
         
    if S in path_SSR: 
        SSR1=path_SSR
        SSR2=exp(-M_SSR.d[t].value)
        D_adj=adjusted_dictionary(M_SSR,D_SSR,R,S,T)
        SSR3=exp(-path_length(M_SSR,path_SSR,D_adj,S,T))
        #SSR3=exp(-(M_SSR.Obj.expr())) #THIS IS PROBABLY WRONG if there are multiple paths
    else: 
        SSR1="N/A"
        SSR2= "N/A"
        SSR3="N/A"
    path_MST=return_path_asymmetric_multiple_ST(M_MST,D,N,S,T)
    path_MSTR=return_path_asymmetric_multiple_ST(M_MSTR,D,N,S,T)
    '''
    evad=0
    for (i,j) in D.keys():
        (l,r)=D[(i,j)]
        evad=evad + l*M_MSTR.STBlocks[(S,T)].z[(i,j)].value + (l+r)*M_MSTR.STBlocks[(S,T)].w[(i,j)].value
    
    evad=evad - sum(M_MSTR.STBlocks[(S,T)].t[i].value*M_MSTR.STBlocks[(S,T)].t[i].value for i in M_MSTR.V)
    '''
    #Fix these two to work with M.STBlocks[(S,T)].t etc etc 
    D_adj=adjusted_dictionary(M_MSTR,D,R,S,T)
    evad=path_length(M_MSTR,path_MSTR,D_adj,S,T)
    
    table.add_row([(S,T), SS1, SS2, SSR1, SSR2, SSR3,path_MST,exp(-M_MST.STBlocks[(S,T)].d[T].value), path_MSTR,exp(-M_MSTR.STBlocks[(S,T)].d[T].value), exp(-evad)])


print(f'Budget={B}')
print(table)

