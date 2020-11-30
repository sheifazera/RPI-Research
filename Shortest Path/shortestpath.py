# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 14:20:34 2020

@author: spunlag
"""
from pyomo.environ import *
import numpy as np
from numpy import linalg as LA
import random
infty=float('inf')
#s is the starting node
#t is the ending node
#N is a set of nodes including s and t
#D is a dictionary mapping (u,v) edges to the length of that edge
#if (u,v) is in the dictionary, do not include the edge (v,u)

def create_shortest_path(D,N,s,t):  #D must provide all edges (so undirected edges should be put twice but with nodes reversed)
    M=ConcreteModel()
    #Check if s and t are in the set of nodes
    if s not in N:
        raise Exception('{s} is not in the provided set of nodes')
    if t not in N:
        raise Exception('{t} is not in the provided set of nodes')
    
    M.d=Var(N,within=NonNegativeReals) #distance from s to a node
    M.Obj=Objective(expr=M.d[t], sense=maximize)
    M.c0=Constraint(expr= M.d[s]==0) #distance from s to s is 0
    
    def C(M,i,j):
        return M.d[j] <= M.d[i] + D[(i,j)]
    M.c=Constraint(D.keys(),rule=C) #Each directed edge should have a constraint           
    M.dual=Suffix(direction=Suffix.IMPORT) #Request the dual variables back from the solver
    return M

def return_path(M,D,N,s,t):
    N_in_path=[s]
    n=s
    nold=infty
    while (n != t):
        if n==nold: #If we went through the loop without finding an adjacent node
            raise Exception('False path: no adjacent node has positive dual value!')
        for j in N-set(N_in_path): 
            if (n,j) in D.keys():
                if M.dual[M.c[(n,j)]] >0:
                    N_in_path.append(j)
                    nold=n
                    n=j
                    break   
    return N_in_path

def create_shortest_path_interdiction(D,N,s,t, B): #D={(u,v):(l,r)}
    M=ConcreteModel()
    #Check if s and t are in the set of nodes
    if s not in N:
        raise Exception('{s} is not in the provided set of nodes')
    if t not in N:
        raise Exception('{t} is not in the provided set of nodes') 
    M.d=Var(N,within=NonNegativeReals) #distance from s to a node
    M.x=Var(D.keys(),within=Binary) #whether or not an edge is interdicted
    
    M.Obj=Objective(expr=M.d[t], sense=maximize)
    M.c0=Constraint(expr= M.d[s]==0) #distance from s to s is 0
    M.budget=Constraint(expr = sum(M.x[(u,v)] for (u,v) in D.keys()) <= B)
    def C(M,i,j):
        (l,r)=D[(i,j)]
        return M.d[j] <= M.d[i] + l + r*M.x[(i,j)]
    M.c=Constraint(D.keys(),rule=C) #Each directed edge should have a constraint  
         
    M.dual=Suffix(direction=Suffix.IMPORT) #Request the dual variables back from the solver
    return M

def create_shortest_path_robust_ellipsoid(D,N,P,s,t,udim): #P is the weight matrix  P={(u,v,j):....}
    ''' max d_t
        s.t. 
        d_j \leq d_i + (l_Pu)_ij for all edges (i,j) in E
        ||u|| <= 1
        
        which is equivalent to:
        
        max d_t
        s.t. 
        d_j \leq d_i + l_ij + t_ij
        t=Pu
        u^Tu \leq w^2 (cone)
        w==1 
    '''
    M=ConcreteModel()
    if len(D.keys())*udim != len(P):
        raise Exception('P is not appropriately dimensioned')
    #Check if s and t are in the set of nodes
    if s not in N:
        raise Exception('{s} is not in the provided set of nodes')
    if t not in N:
        raise Exception('{t} is not in the provided set of nodes') 
    M.d=Var(N,within=PositiveReals) #distance from s to a node
    
    M.U=RangeSet(1,udim)
    
    #M.P=Param(D.keys(),M.U, initialize=P)
    M.u=Var(M.U,within=Reals)
    M.t=Var(D.keys(),within=Reals)
    M.Obj=Objective(expr=M.d[t], sense=maximize)
    M.c0=Constraint(expr= M.d[s]==0) #distance from s to s is 0
    
    def c1(M,u,v):
        value=sum(P[u,v,j-1]*M.u[j] for j in M.U) 
        return  value == M.t[(u,v)]
    M.c1=Constraint(D.keys(),rule=c1)
    
    def C(M,i,j):
        return M.d[j] <= M.d[i] + D[(i,j)] + M.t[(i,j)]
    M.c=Constraint(D.keys(),rule=C)
    
    M.w=Var(within=NonNegativeReals)
    def c2(M):
        value=sum(M.u[i]*M.u[i] for i in M.U)
        return value <= M.w*M.w
    M.c2=Constraint(rule=c2)
    M.c3=Constraint(expr=M.w==1)
    
    M.dual=Suffix(direction=Suffix.IMPORT) #Request the dual variables back from the solver
    return M

def create_shortest_path_robust_ellipsoid_certain_interdiction(D,N,P,s,t,udim, B): #P is the weight matrix  P={(u,v,j):....} D=(i,j):(l,r)
    M=ConcreteModel()
    if len(D.keys())*udim != len(P):
        raise Exception('P is not appropriately dimensioned')
    #Check if s and t are in the set of nodes
    if s not in N:
        raise Exception('{s} is not in the provided set of nodes')
    if t not in N:
        raise Exception('{t} is not in the provided set of nodes') 
    M.d=Var(N,within=PositiveReals) #distance from s to a node
    M.x=Var(D.keys(),within=Binary) #whether or not an edge is interdicted
    M.budget=Constraint(expr = sum(M.x[(u,v)] for (u,v) in D.keys()) <= B)
    M.U=RangeSet(1,udim)
    #M.P=Param(D.keys(),M.U, initialize=P)
    M.u=Var(M.U,within=Reals)
    M.t=Var(D.keys(),within=Reals)
    M.Obj=Objective(expr=M.d[t], sense=maximize)
    M.c0=Constraint(expr= M.d[s]==0) #distance from s to s is 0
    
    def c1(M,u,v):
        value=sum(P[u,v,j-1]*M.u[j] for j in M.U) 
        return  value == M.t[(u,v)]
    M.c1=Constraint(D.keys(),rule=c1)
    
    def C(M,i,j):
        (l,r)=D[(i,j)]
        return M.d[j] <= M.d[i] + l + M.t[(i,j)] + r*M.x[(i,j)]
    M.c=Constraint(D.keys(),rule=C)
    
    M.w=Var(within=NonNegativeReals)
    def c2(M):
        value=sum(M.u[i]*M.u[i] for i in M.U)
        return value <= M.w*M.w
    M.c2=Constraint(rule=c2)
    M.c3=Constraint(expr=M.w==1)
    
    M.dual=Suffix(direction=Suffix.IMPORT) #Request the dual variables back from the solver
    return M












def return_path_robust(M,D,N,P,s,t,udim):
    N_in_path=[s]
    n=s
    nold=infty
    tol=1/len(N)
    while (n != t):
        if n==nold: #If we went through the loop without finding an adjacent node
            raise Exception('False path: no adjacent node has positive dual value!')
        for j in N-set(N_in_path): 
            if (n,j) in D.keys():
                if M.dual[M.c[(n,j)]] >tol:
                    N_in_path.append(j)
                    nold=n
                    n=j
                    break   
    return N_in_path


def return_path_interdiction(M,D,N,s,t,opt):
    D_interdicted={}
    for (i,j) in D.keys():
        (l,r)=D[(i,j)]
        if M.x[(i,j)]==1:
            D_interdicted[(i,j)]=l+r
        else:
            D_interdicted[(i,j)]=l
    M_interdicted=create_shortest_path(D_interdicted,N,s,t)
    opt.solve(M_interdicted)
    path=return_path(M_interdicted,D_interdicted,N,s,t)
    return path

def return_path_interdiction_robust(M,D,N,P,s,t,udim,opt): 

    D_interdicted={}
    for (i,j) in D.keys():
        (l,r)=D[(i,j)]
        if M.x[(i,j)]==1:
            D_interdicted[(i,j)]=l + r + sum(P[(i,j,k)]*M.u[k+1].value for k in range(0,udim))
        else:
            D_interdicted[(i,j)]=l +     sum(P[(i,j,k)]*M.u[k+1].value for k in range(0,udim))
    M_interdicted=create_shortest_path(D_interdicted,N,s,t)
    opt.solve(M_interdicted)
    path=return_path(M_interdicted,D_interdicted,N,s,t)
    return path


def create_shortest_path_robust_ellipsoid_uncertain_interdiction(D,N,P,R, s,t,udim,vdim, B):
    M=ConcreteModel()
    if len(D.keys())*udim != len(P):
        raise Exception('P is not appropriately dimensioned')
    #Check if s and t are in the set of nodes
    if s not in N:
        raise Exception('{s} is not in the provided set of nodes')
    if t not in N:
        raise Exception('{t} is not in the provided set of nodes') 
    M.d=Var(N,within=PositiveReals) #distance from s to a node
    M.x=Var(D.keys(),within=Binary) #whether or not an edge is interdicted
    M.p=Var(D.keys(),within=NonNegativeReals) #p_ij=x_ij*s_ij
    M.s=Var(D.keys(),within=Reals) #s=Rv
    M.budget=Constraint(expr = sum(M.x[(u,v)] for (u,v) in D.keys()) <= B)
    M.U=RangeSet(1,udim)
    M.V=RangeSet(1,vdim)
    #M.P=Param(D.keys(),M.U, initialize=P)
    M.u=Var(M.U,within=Reals)
    M.v=Var(M.V,within=Reals)
    M.t=Var(D.keys(),within=Reals)
    M.Obj=Objective(expr=M.d[t], sense=maximize)
    M.c0=Constraint(expr= M.d[s]==0) #distance from s to s is 0
    
    def c1(M,u,v):
        value=sum(P[u,v,j-1]*M.u[j] for j in M.U) 
        return  value == M.t[(u,v)]
    M.c1=Constraint(D.keys(),rule=c1)
    
    def C(M,i,j):
        (l,r)=D[(i,j)]
        return M.d[j] <= M.d[i] + l + M.t[(i,j)] + r*M.x[(i,j)] + M.p[(i,j)]
    M.c=Constraint(D.keys(),rule=C)
    
    M.w=Var(within=NonNegativeReals)
    def c2(M):
        value=sum(M.u[i]*M.u[i] for i in M.U)
        return value <= M.w*M.w
    M.c2=Constraint(rule=c2)
    M.c3=Constraint(expr=M.w==1)
    def c4(M,i,j):
        value=sum(R[i,j,k-1]*M.v[k] for k in M.V) 
        return  value == M.s[(i,j)] #Rv=s
    M.c4=Constraint(D.keys(),rule=c4)
    def c5(M):
        value=sum(M.v[i]*M.v[i] for i in M.V)
        return value <= M.w*M.w #v^Tv\leq w^2, w==1
    M.c5=Constraint(rule=c5)
    
    R_array=np.empty([len(D.keys()),vdim],dtype=float)
    l=0
    for (i,j) in D.keys():
        for k in M.V:
            R_array[l,k-1]=R[(i,j,k-1)]
        l=l+1    
    R_norm=LA.norm(R_array,2)
    
    def c6(M,i,j):
        value=R_norm*M.x[(i,j)]+ M.s[(i,j)] - R_norm
        return M.p[(i,j)] >= value
    M.c6=Constraint(D.keys(), rule=c6)
    def c7(M,i,j):
        value=R_norm*M.x[(i,j)]
        return M.p[(i,j)] <= value
    M.c7=Constraint(D.keys(),rule=c7)
    def c8(M,i,j):
        return M.p[(i,j)] <= M.s[(i,j)]
    M.c8=Constraint(D.keys(),rule=c8)
    
    M.dual=Suffix(direction=Suffix.IMPORT) #Request the dual variables back from the solver
    
    return M


def return_path_interdiction_robust_uncertain(M,D,N,P,R,s,t,udim,vdim,opt): 
    D_interdicted={}
    for (i,j) in D.keys():
        (l,r)=D[(i,j)]
        if M.x[(i,j)]==1:
            D_interdicted[(i,j)]=l + r + sum(P[(i,j,k)]*M.u[k+1].value for k in range(0,udim)) + sum(R[(i,j,k)]*M.v[k+1].value for k in range(0,vdim))
        else:
            D_interdicted[(i,j)]=l +     sum(P[(i,j,k)]*M.u[k+1].value for k in range(0,udim))
    M_interdicted=create_shortest_path(D_interdicted,N,s,t)
    opt.solve(M_interdicted)
    path=return_path(M_interdicted,D_interdicted,N,s,t)
    return path


def create_asymmetric_uncertainty_shortest_path_interdiction(D,N,R,s,t,vdim,B):
    M=ConcreteModel()
    if len(D.keys())*vdim != len(R):
        raise Exception('R is not appropriately dimensioned')
    #Check if s and t are in the set of nodes
    if s not in N:
        raise Exception('{s} is not in the provided set of nodes')
    if t not in N:
        raise Exception('{t} is not in the provided set of nodes') 
    M.d=Var(N,within=PositiveReals) #distance from s to a node
    M.x=Var(D.keys(),within=Binary) #whether or not an edge is interdicted
    M.budget=Constraint(expr = sum(M.x[(u,v)] for (u,v) in D.keys()) <= B)
    M.V=RangeSet(1,vdim)
    M.v=Var(M.V,within=Reals)
    M.t=Var(M.V,within=Reals)
    M.w=Var(D.keys(),within=NonNegativeReals)
    M.z=Var(D.keys(),within=NonNegativeReals)
    M.s=Var(within=NonNegativeReals)
    
    def Obj(M):
        value=0
        for (i,j) in D.keys():
            (l,r)=D[(i,j)]
            value=value + l*M.z[(i,j)] + (l+r)*M.w[(i,j)]
        value=value -M.s
        return value
    M.Obj=Objective(rule=Obj,sense=maximize)
    
    
    def C(M,i,j):
        (l,r)=D[(i,j)]
        return M.d[j] <= M.d[i] + l + r*M.x[(i,j)]
    M.c=Constraint(D.keys(),rule=C) #Distance-based shortest path
    M.c0=Constraint(expr = M.d[s]==0) #Distance-based shortest path
    def c1(M,k):
        value=sum(R[(i,j,k-1)]*M.w[(i,j)] for (i,j) in D.keys())
        return M.t[k]==value
    M.c1=Constraint(M.V, rule=c1) # R^Tw=t
    
    M.c2=Constraint(expr= sum(M.t[i]*M.t[i] for i in M.V)<= M.s*M.s) # t^T t <= s^2
    def c3(M,i,j):
        return M.z[(i,j)] + M.x[(i,j)] <=1 
    M.c3=Constraint(D.keys(),rule=c3) #zk + xk <= 1
    def c4(M,i,j):
        return M.w[(i,j)] - M.x[(i,j)] <= 0
    M.c4=Constraint(D.keys(),rule=c4) #wk - xk <= 0
    def c5(M):
        value=0
        for (i,j) in D.keys():
            (l,r)=D[(i,j)]
            value=value+ l*M.z[(i,j)] + (l+r)*M.w[(i,j)]
        return M.d[t]==value
    M.c5=Constraint(rule=c5)
    
    def c6(M,J):
        if J==s:
            rhs=1
        elif J==t:
            rhs=-1
        else:
            rhs=0
        RS=0
        FS=0
        for (i,j) in D.keys():
            if j==J:
                RS= RS + M.w[(i,j)] + M.z[(i,j)]
            elif i==J:
                FS= FS + M.w[(i,j)] + M.z[(i,j)]
        return FS-RS == rhs 
    M.c6=Constraint(N,rule=c6)
  
    return M


def return_path_asymmetric(M,D,N,s,t):
    tol=1/max(N)
    N_in_path=[s]
    n=s
    nold=infty
    while (n != t):
        if n==nold: #If we went through the loop without finding an adjacent node
            raise Exception('False path: no adjacent node has positive dual value!')
        for j in N-set(N_in_path): 
            if (n,j) in D.keys():
                if M.z[(n,j)].value >tol or M.w[(n,j)].value >tol:
                    N_in_path.append(j)
                    nold=n
                    n=j
                    break   
    return N_in_path





#WIP
 
def create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST(D,N,R,I,vdim,B):
    M=ConcreteModel()
    if len(D.keys())*vdim != len(R):
        raise Exception('R is not appropriately dimensioned')   
    M.x=Var(D.keys(),within=Binary) #whether or not an edge is interdicted
    M.budget=Constraint(expr = sum(M.x[(i,j)] for (i,j) in D.keys()) <= B)
    M.V=RangeSet(1,vdim)
    M.m=Var(within=Reals)
    
    M.Obj=Objective(expr=M.m, sense=maximize)
    
    M.STBlocks=Block(I)
    
    #Loop for all possible (S,T) pairs
    for (S,T) in I:  
        M.STBlocks[(S,T)].d=Var(N,within=PositiveReals)
        M.STBlocks[(S,T)].t=Var(M.V,within=Reals)
        M.STBlocks[(S,T)].w=Var(D.keys(),within=NonNegativeReals)
        M.STBlocks[(S,T)].z=Var(D.keys(),within=NonNegativeReals)
        M.STBlocks[(S,T)].s=Var(within=NonNegativeReals)
    #Constraint Rules Here
        def C(M,i,j):
            (l,r)=D[(i,j)]
            return M.d[j] <= M.d[i] + l + r*M.x[(i,j)]
        def c1(M,k):
            value=sum(R[(i,j,k-1)]*M.w[(S,T),(i,j)] for (i,j) in D.keys())
            return M.t[S,T,k]==value
        def c3(M,i,j):
            return M.z[(S,T),(i,j)] + M.x[(i,j)] <=1 
        def c4(M,i,j):
            return M.w[(S,T),(i,j)] - M.x[(i,j)] <= 0
        def c5(M):
            value=0
            for (i,j) in D.keys():
                (l,r)=D[(i,j)]
                value=value+ l*M.z[(S,T),(i,j)] + (l+r)*M.w[(S,T),(i,j)]
            return M.d[(S,T),T]==value
        def c6(M,J):
            if J==S:
                rhs=1
            elif J==T:
                rhs=-1
            else:
                rhs=0
            RS=0
            FS=0
            for (i,j) in D.keys():
                if j==J:
                    RS= RS + M.w[(S,T),(i,j)] + M.z[(S,T),(i,j)]
                elif i==J:
                    FS= FS + M.w[(S,T),(i,j)] + M.z[(S,T),(i,j)]
            return FS-RS == rhs 
        
        def c7(M):
            value=0
            for (i,j) in D.keys():
                (l,r)=D[(i,j)]
                value=value + l*M.z[(S,T),(i,j)] + (l+r)*M.w[(S,T),(i,j)]
            value=value -M.s[(S,T)]
            
            return value >= M.m
    

        M.STBlocks[(S,T)].c=Constraint(D.keys(),rule=C) #Distance-based shortest path
        M.STBlocks[(S,T)].c0=Constraint(expr = M.d[(S,T),S]==0) #Distance-based shortest path
        M.STBlocks[(S,T)].c1=Constraint(M.V, rule=c1) # R^Tw=t
        M.STBlocks[(S,T)].c2=Constraint(expr= sum(M.t[(S,T),i]*M.t[(S,T),i] for i in M.V)<= M.s[(S,T)]*M.s[(S,T)]) # t^T t <= s^2
        M.STBlocks[(S,T)].c3=Constraint(D.keys(),rule=c3) #zk + xk <= 1
        M.STBlocks[(S,T)].c4=Constraint(D.keys(),rule=c4) #wk - xk <= 0
        M.STBlocks[(S,T)].c5=Constraint(rule=c5)
        M.STBlocks[(S,T)].c6=Constraint(N,rule=c6)
        M.STBlocks[(S,T)].c7=Constraint(rule=c7)
 
    return M














'''
def return_path(M,D,N,s,t):
    N_in_path=[t] #Use lists not sets to maintain order (going backwards)
    N_not_in_path=[]
    n=t 
    while (n != s): 
        for u in N: #Do NOT consider nodes already visited or nodes that we've checked cannot reach the end properly
            if u not in set(N_not_in_path).union(set(N_in_path)):
                
                if (u,n) in D.keys():
                    if M.d[u].value==M.d[n].value-D[(u,n)]:
                        N_in_path.append(u)
                        n=u
                        break
                elif u==max(N): #if I didn't break and I reach the end of the nodes, this node shouldn't be in the path
                    N_in_path.remove(n) #remove most recent
                    N_not_in_path.append(n)
                    n=N_in_path[-1] #change the new n to the previous node in the current path
                
    return N_in_path

def return_path_robust(M,D,N,P,s,t,udim):
    l={}
    for (i,j) in D.keys():
        l[(i,j)]=D[(i,j)] + sum(P[(i,j,k)]*M.u[k+1].value for k in range(0,udim))
    N_in_path=[t] #Use lists not sets to maintain order (going backwards)
    N_not_in_path=[]
    n=t 
    while (n != s): 
        for u in N: #Do NOT consider nodes already visited or nodes that we've checked cannot reach the end properly
            if u not in set(N_not_in_path).union(set(N_in_path)):
                
                if (u,n) in D.keys():
                    if M.d[u].value==M.d[n].value-l[(u,n)]:
                        N_in_path.append(u)
                        n=u
                        break
                elif u==max(N): #if I didn't break and I reach the end of the nodes, this node shouldn't be in the path
                    N_in_path.remove(n) #remove most recent
                    N_not_in_path.append(n)
                    n=N_in_path[-1] #change the new n to the previous node in the current path
                
    return N_in_path
'''


'''
This is deprecated: create_shortest_path now requires directed arcs
def create_directed_shortest_path(D,N,s,t): #(u,v) goes FROM U TO V, (v,u) goes FROM V TO U and should both be included if the arc goes both ways
    M=ConcreteModel()
    #Check if s and t are in the set of nodes
    if s not in N:
        raise Exception('{s} is not in the provided set of nodes')
    if t not in N:
        raise Exception('{t} is not in the provided set of nodes')
    
    M.d=Var(N,within=PositiveReals) #distance from s to a node
    M.Obj=Objective(expr=M.d[t], sense=maximize)
    M.c0=Constraint(expr= M.d[s]==0) #distance from s to s is 0
    M.c=ConstraintList()
    for node in N:
        for _Eid, L in D.items():
            (u,v)=_Eid
            if v==node:
                M.c.add(M.d[v] <= M.d[u] +L ) 
    return M
'''


