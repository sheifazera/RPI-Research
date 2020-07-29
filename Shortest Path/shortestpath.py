# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 14:20:34 2020

@author: spunlag
"""
from pyomo.environ import *
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