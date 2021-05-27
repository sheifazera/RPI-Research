# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 10:20:35 2021

@author: sheifazera
"""

from pyomo.environ import *
import numpy as np
from numpy import linalg as LA
import random
import networkx as nx
infty=float('inf')
from itertools import islice



#Shortest Path
def create_shortest_path_nx(G,s,t):  #G is a directed NetworkX graph
    M=ConcreteModel()
    #Check if s and t are in the set of nodes
    if s not in list(G.nodes):
        raise Exception('{s} is not in the provided set of nodes')
    if t not in list(G.nodes):
        raise Exception('{t} is not in the provided set of nodes')
    
    M.d=Var(set(G.nodes),within=NonNegativeReals) #distance from s to a node
    M.Obj=Objective(expr=M.d[t], sense=maximize)
    M.c0=Constraint(expr= M.d[s]==0) #distance from s to s is 0
    
    def C(M,i,j):
        
        return M.d[j] <= M.d[i] + G[i][j]['length']
    M.c=Constraint(set(G.edges),rule=C) #Each directed edge should have a constraint           
    M.dual=Suffix(direction=Suffix.IMPORT) #Request the dual variables back from the solver
    return M



#Shortest Path with Interdiction
def create_shortest_path_interdiction_nx(G,s,t, B): #G is a directed NetworkX graph with length and interdicted_length
    M=ConcreteModel()
    #Check if s and t are in the set of nodes
    if s not in list(G.nodes):
        raise Exception('{s} is not in the provided set of nodes')
    if t not in list(G.nodes):
        raise Exception('{t} is not in the provided set of nodes') 
    M.d=Var(set(G.nodes),within=NonNegativeReals) #distance from s to a node
    M.x=Var(set(G.edges),within=Binary) #whether or not an edge is interdicted
    
    M.Obj=Objective(expr=M.d[t], sense=maximize)
    M.c0=Constraint(expr= M.d[s]==0) #distance from s to s is 0
    M.budget=Constraint(expr = sum(M.x[(u,v)] for (u,v) in set(G.edges)) <= B)
    def C(M,i,j):
        return M.d[j] <= M.d[i] + G[i][j]['length'] + M.x[(i,j)]* G[i][j]['interdicted_length']
    M.c=Constraint(set(G.edges),rule=C) #Each directed edge should have a constraint       
    return M




#Shortest Path with Asymmetric Robust Uncertainty 
def create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R,s,t,vdim,B):
    M=ConcreteModel()
    if len(set(G.edges))*vdim != len(R):
        raise Exception('R is not appropriately dimensioned')
    #Check if s and t are in the set of nodes
    if s not in list(G.nodes):
        raise Exception('{s} is not in the provided set of nodes')
    if t not in list(G.nodes):
        raise Exception('{t} is not in the provided set of nodes') 
    M.d=Var(set(G.nodes),within=PositiveReals) #distance from s to a node
    M.x=Var(set(G.edges),within=Binary) #whether or not an edge is interdicted
    M.budget=Constraint(expr = sum(M.x[(u,v)] for (u,v) in set(G.edges)) <= B)
    M.V=RangeSet(1,vdim)
    M.v=Var(M.V,within=Reals)
    M.t=Var(M.V,within=Reals)
    M.w=Var(set(G.edges),within=NonNegativeReals)
    M.z=Var(set(G.edges),within=NonNegativeReals)
    M.s=Var(within=NonNegativeReals)
    
    def Obj(M):
        value=0
        for (i,j) in set(G.edges):
            l=G[i][j]['length']
            r=G[i][j]['interdicted_length']
            value=value + l*M.z[(i,j)] + (l+r)*M.w[(i,j)]
        value=value -M.s
        return value
    M.Obj=Objective(rule=Obj,sense=maximize)
    
    
    def C(M,i,j):
        return M.d[j] <= M.d[i] + G[i][j]['length'] + M.x[(i,j)]* G[i][j]['interdicted_length']
    M.c=Constraint(set(G.edges),rule=C) #Distance-based shortest path
    M.c0=Constraint(expr = M.d[s]==0) #Distance-based shortest path
    def c1(M,k):
        value=sum(R[(i,j,k-1)]*M.w[(i,j)] for (i,j) in set(G.edges))
        return M.t[k]==value
    M.c1=Constraint(M.V, rule=c1) # R^Tw=t
    
    M.c2=Constraint(expr= sum(M.t[i]*M.t[i] for i in M.V)<= M.s*M.s) # t^T t <= s^2
    def c3(M,i,j):
        return M.z[(i,j)] + M.x[(i,j)] <=1 
    M.c3=Constraint(set(G.edges),rule=c3) #zk + xk <= 1
    def c4(M,i,j):
        return M.w[(i,j)] - M.x[(i,j)] <= 0
    M.c4=Constraint(set(G.edges),rule=c4) #wk - xk <= 0
    def c5(M):
        value=0
        for (i,j) in set(G.edges):
            l=G[i][j]['length']
            r=G[i][j]['interdicted_length']
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
        for (i,j) in set(G.out_edges(J)):
            FS= FS + M.w[(i,j)] + M.z[(i,j)]    
        for (i,j) in set(G.in_edges(J)):
            RS= RS + M.w[(i,j)] + M.z[(i,j)]
        return FS-RS == rhs 
    M.c6=Constraint(set(G.nodes),rule=c6)
  
    return M



#Shortest Path with Asymmetric Robust Uncertainty for multiple (S,T) pairs

def create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R,I,vdim,B):
    M=ConcreteModel()
    if  len(set(G.edges))*vdim != len(R):
        raise Exception('R is not appropriately dimensioned')   
    M.x=Var(set(G.edges),within=Binary) #whether or not an edge is interdicted
    M.budget=Constraint(expr = sum(M.x[(u,v)] for (u,v) in set(G.edges)) <= B)
    M.V=RangeSet(1,vdim)
    M.m=Var(within=Reals)
    
    M.Obj=Objective(expr=M.m, sense=maximize)
    
    M.STBlocks=Block(I)
    
    #Loop for all possible (S,T) pairs
    for (S,T) in I:  
        M.STBlocks[(S,T)].d=Var(set(G.nodes),within=PositiveReals)
        M.STBlocks[(S,T)].t=Var(M.V,within=Reals)
        M.STBlocks[(S,T)].w=Var(set(G.edges),within=NonNegativeReals)
        M.STBlocks[(S,T)].z=Var(set(G.edges),within=NonNegativeReals)
        M.STBlocks[(S,T)].s=Var(within=NonNegativeReals)
    #Constraints
        M.STBlocks[(S,T)].C=ConstraintList()
        for (i,j) in set(G.edges):
            l=G[i][j]['length']
            r=G[i][j]['interdicted_length']
            M.STBlocks[(S,T)].C.add(M.STBlocks[(S,T)].d[j] <= M.STBlocks[(S,T)].d[i] + l + r*M.x[(i,j)])
        M.STBlocks[(S,T)].c0=Constraint(expr = M.STBlocks[(S,T)].d[S]==0)
        M.STBlocks[(S,T)].c1=ConstraintList()
        for k in M.V:
            M.STBlocks[(S,T)].c1.add(M.STBlocks[(S,T)].t[k]==sum(R[(i,j,k-1)]*M.STBlocks[(S,T)].w[(i,j)] for (i,j) in set(G.edges)))
        M.STBlocks[(S,T)].c2=Constraint(expr=sum(M.STBlocks[(S,T)].t[i]*M.STBlocks[(S,T)].t[i] for i in M.V)<= M.STBlocks[(S,T)].s*M.STBlocks[(S,T)].s)
        M.STBlocks[(S,T)].c3=ConstraintList()
        for (i,j) in set(G.edges):
            M.STBlocks[(S,T)].c3.add(M.STBlocks[(S,T)].z[(i,j)] + M.x[(i,j)] <=1 )
        M.STBlocks[(S,T)].c4=ConstraintList()
        for (i,j) in set(G.edges):
            M.STBlocks[(S,T)].c4.add(M.STBlocks[(S,T)].w[(i,j)] - M.x[(i,j)] <= 0)
        value=0
        for (i,j) in set(G.edges):
            l=G[i][j]['length']
            r=G[i][j]['interdicted_length']
            value=value+ l*M.STBlocks[(S,T)].z[(i,j)] + (l+r)*M.STBlocks[(S,T)].w[(i,j)]
        M.STBlocks[(S,T)].c5=Constraint(expr=M.STBlocks[(S,T)].d[T]==value)
        M.STBlocks[(S,T)].c6=ConstraintList()
        for J in set(G.nodes):
            if J==S:
                rhs=1
            elif J==T:
                rhs=-1
            else:
                rhs=0
            RS=0
            FS=0
            for (i,j) in set(G.out_edges(J)):
                FS= FS + M.STBlocks[(S,T)].w[(i,j)] + M.STBlocks[(S,T)].z[(i,j)]
            for (i,j) in set(G.in_edges(J)):
                RS= RS + M.STBlocks[(S,T)].w[(i,j)] + M.STBlocks[(S,T)].z[(i,j)]
                   
            M.STBlocks[(S,T)].c6.add(FS-RS == rhs )        
        
    
        value=0
        for (i,j) in set(G.edges):
            l=G[i][j]['length']
            r=G[i][j]['interdicted_length']
            value=value + l*M.STBlocks[(S,T)].z[(i,j)] + (l+r)*M.STBlocks[(S,T)].w[(i,j)]
        value=value -M.STBlocks[(S,T)].s
        M.STBlocks[(S,T)].c7=Constraint(expr=value >= M.m)
 
    return M


#Return path




def return_paths(M,G,R,S,T):
    tol=1e-5
    G_adj=nx.DiGraph()
    norm=np.sqrt(sum(M.t[i].value*M.t[i].value for i in M.V))
    for (i,j) in set(G.edges):
        if (M.w[(i,j)].value >= 1/len(list(G.nodes))) or (M.z[(i,j)].value >= 1/len(list(G.nodes))):
            l=G[i][j]['length']
            r=G[i][j]['interdicted_length']
            if M.x[(i,j)].value >=0.9:
                length=l+r
                if norm>tol: # Tune this cut off value
                    adjustment=(-1/norm)*(sum(R[(i,j,k-1)]*M.t[k].value for k in M.V))
                    length=length+adjustment
            else:
                length=l
            G_adj.add_edge(i,j,adjusted_length=length)
    #Finished creating a network with only edges that are traversed in the problem with adjusted edges lengths
    #Use Network_X to find ALL shortest paths
    paths=nx.all_shortest_paths(G_adj, source=S, target=T, weight='adjusted_length')
    paths=list(paths)
    lengths=nx.shortest_path_length(G_adj, source=S, target=T, weight='adjusted_length', method='dijkstra')
    return paths,lengths

def return_interdicted_arcs(M,G):
    interdicted=[]
    for (i,j) in set(G.edges):
        if M.x[(i,j)].value >= 0.9:
            interdicted.append((i,j))
    return interdicted
        
def return_paths_multiple_ST(M,G,R,S,T):
    tol=1e-5
    G_adj=nx.DiGraph()
    norm=np.sqrt(sum(M.STBlocks[(S,T)].t[i].value*M.STBlocks[(S,T)].t[i].value for i in M.V))
    for (i,j) in set(G.edges):
        if (M.STBlocks[(S,T)].w[(i,j)].value >= 1/len(list(G.nodes))) or (M.STBlocks[(S,T)].z[(i,j)].value >= 1/len(list(G.nodes))):
            l=G[i][j]['length']
            r=G[i][j]['interdicted_length']
            if M.x[(i,j)].value >=0.9:
                length=l+r
                if norm>tol: # Tune this cut off value
                    adjustment=(-1/norm)*(sum(R[(i,j,k-1)]*M.STBlocks[(S,T)].t[k].value for k in M.V))
                    length=length+adjustment
            else:
                length=l
            G_adj.add_edge(i,j,adjusted_length=length)
    #Finished creating a network with only edges that are traversed in the problem with adjusted edges lengths
    #Use Network_X to find ALL shortest paths
    paths=nx.all_shortest_paths(G_adj, source=S, target=T, weight='adjusted_length')
    paths=list(paths)
    lengths=nx.shortest_path_length(G_adj, source=S, target=T, weight='adjusted_length', method='dijkstra')
    return paths, lengths


def similar_length(M,G,R,S,T,cutoff): #cutoff=percentage of shortest path longer ie 110% would be cutoff=1.10
    tol=1e-5    
    (paths, lengths)=return_paths(M,G,R,S,T)   
    
    G_adj=nx.Graph()
    norm=np.sqrt(sum(M.t[i].value*M.t[i].value for i in M.V))
    for (i,j) in G.edges:
        l=G[i][j]['length']
        r=G[i][j]['interdicted_length']
        if M.x[(i,j)].value >=0.9:
            length=l+r
            if norm>tol:
                adjustment=(-1/norm)*(sum(R[(i,j,k-1)]*M.t[k].value for k in M.V))
                length=length+adjustment
        else:
            length=l
        if G_adj.has_edge(j,i):
            #Check which length to use
            if G_adj[j][i]['adjusted_length'] <= length: #Shorter Lengths indicate that the edge has been interdicted so we shoudl use that value for both directions
                G_adj[i][j]['adjusted_length']=length
        else:
            G_adj.add_edge(i,j,adjusted_length=length)
    
    similar_paths=[]
    all_paths=list(islice(nx.shortest_simple_paths(G_adj, S, T, weight='adjusted_length'), 30))
    
    for path in all_paths:
        if nx.path_weight(G_adj, path, weight='adjusted_length') <= cutoff*lengths:
            similar_paths.append(path)
        else:
            break
    no_of_similar_paths=len(similar_paths)
    return (similar_paths, no_of_similar_paths)
    
    