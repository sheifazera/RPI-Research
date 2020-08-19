# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 14:20:34 2020

@author: spunlag
"""
from pyomo.environ import *
import numpy as np
import random
#s is the starting node
#t is the ending node
#N is a set of nodes including s and t
#D is a dictionary mapping (u,v) edges to the length of that edge
#if (u,v) is in the dictionary, do not include the edge (v,u)

def create_shortest_path(D,N,s,t): #This is just shortest path no interdiction
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
            if u==node:
                M.c.add(M.d[u] <= M.d[v] + L)
            if v==node:
                M.c.add(M.d[v] <= M.d[u] +L ) 
    return M


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

def create_shortest_path_interdiction(D,N,s,t, B): #D={(u,v):(l,r)}
    M=ConcreteModel()
    #Check if s and t are in the set of nodes
    if s not in N:
        raise Exception('{s} is not in the provided set of nodes')
    if t not in N:
        raise Exception('{t} is not in the provided set of nodes') 
    M.d=Var(N,within=PositiveReals) #distance from s to a node
    M.x=Var(D.keys(),within=Binary) #whether or not an edge is interdicted
    
    M.Obj=Objective(expr=M.d[t], sense=maximize)
    M.c0=Constraint(expr= M.d[s]==0) #distance from s to s is 0
    M.budget=Constraint(expr = sum(M.x[(u,v)] for (u,v) in D.keys()) <= B)
    M.c=ConstraintList()
    for node in N:
        for _Eid, _param in D.items():
            (u,v)=_Eid
            (L,R)=_param
            if u==node:
                M.c.add(M.d[u] <= M.d[v] + L + R* M.x[(u,v)])
            if v==node:
                M.c.add(M.d[v] <= M.d[u] + L  + R*M.x[(u,v)]) 
    return M

def create_shortest_path_robust_ellipsoid(D,N,P,s,t,udim): #P is the weight matrix  P={(u,v,j):....}
    M=ConcreteModel()
    '''
    if m!=len(N):
        raise Exception('P is not appropriately dimensioned')
    '''
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
        value=sum(P[u,v,j-1]*M.u[j] for j in M.U) #ERROR P is not indexed by the edges, but a range set...
        return  value == M.t[(u,v)]
    M.c1=Constraint(D.keys(),rule=c1)
    
    M.c=ConstraintList()
    for node in N:
        for _Eid, L in D.items():
            (u,v)=_Eid
            if u==node:
                M.c.add(M.d[u] <= M.d[v] + L + M.t[_Eid])
            if v==node:
                M.c.add(M.d[v] <= M.d[u] + L + M.t[_Eid]) 
    M.w=Var(within=NonNegativeReals)
    def c2(M):
        value=sum(M.u[i]*M.u[i] for i in M.U)
        return value <= M.w*M.w
    M.c2=Constraint(rule=c2)
    M.c3=Constraint(expr=M.w==1)
    return M






D={(0,1):22,(0,2):11,(1,3):69,(2,3):120}
N={0,1,2,3}
s=0
t=3

M=create_shortest_path(D,N,s,t)
opt=SolverFactory('gurobi')
results=opt.solve(M)
M.pprint()


D={(0,1):(22,29),(0,2):(11,40),(1,3):(69,22),(2,3):(120,109)}
N={0,1,2,3}
s=0
t=3
B=1
M=create_shortest_path_interdiction(D, N, s, t, B)
results=opt.solve(M)
M.pprint()

#assume u is 3-dimensional so P is 4x3 dimensional also the edges go in the same order they do in the dictionary D...OR ELSE

D={(0,1):22,(0,2):11,(1,3):69,(2,3):120}
udim=3
P_array=np.random.uniform(low=0, high=1, size=(4,udim))
P={}
j=0
for (u,v) in D.keys():
    for i in range(0,udim):
        P[(u,v,i)]=P_array[j,i]
    j=j+1
    
#P={(0,1,0):random.randrange(0,1),(0,1,1):random.randrange(0,1),(0,1,2):random.randrange(0,1),(0,1,3):random.randrange(0,1)}

M_ellipsoid=create_shortest_path_robust_ellipsoid(D, N, P, s, t, udim)