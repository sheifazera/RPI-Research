# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 12:34:08 2020

@author: sheifazera
"""

import time

from pyomo.environ import *
from pao.bilevel import *

from NuclearSmugglingParameters import E, FS, RS, n, e, b, c, p, q, P, ns, nt, W

#A incidence matrix
#n number of nodes
#e number of edges
#b budget
#c cost of edge sensor placement
#ns number of starting nodes
#nt number of ending nodes


m=ConcreteModel()
m.nodeset=RangeSet(1,n)
m.edgeset=RangeSet(1,e)
m.Edgeset=Set(initialize=E)
m.Omega=RangeSet(1,nt*ns)
m.OmegaSet=Set(initialize=W)
m.Node=Block(m.nodeset)

for j in m.nodeset:
    m.Node[j].FS=Set(initialize=FS[j])
    m.Node[j].RS=Set(initialize=RS[j])
    
m.b=Param(initialize=b)
m.c=Param(m.Edgeset,initialize=c)
m.P=Param(m.OmegaSet,initialize=P)
m.p=Param(m.Edgeset,initialize=p)
m.q=Param(m.Edgeset,initialize=q)
m.x=Var(m.Edgeset,within=Binary)
m.h=Var(m.OmegaSet,within=NonNegativeReals)
m.y=Var(m.OmegaSet,m.Edgeset,within=NonNegativeReals)
m.z=Var(m.OmegaSet,m.Edgeset,within=NonNegativeReals)

def Obj(m):
    value=sum(m.P[(i,j)]*m.h[(i,j)] for (i,j) in m.OmegaSet) #h[i,j] is the objective of subproblem (i,j)
    return value
m.Obj=Objective(rule=Obj,sense=minimize)

def C(m):
    value=sum(m.c[(i,j)]*m.x[(i,j)] for (i,j) in m.Edgeset)
    return value <= b
m.C=Constraint(rule=C)

#m.pprint()

m.Sub=Block(m.OmegaSet)
for (s,t) in m.OmegaSet:
    m.Sub[(s,t)].sub=SubModel(fixed=m.x)
    m.Sub[(s,t)].sub.stset=Set(initialize=[s,t])
    def obj(sub):
        value=m.h[(s,t)]
        return value
    m.Sub[(s,t)].sub.obj=Objective(rule=obj, sense=maximize)
    
    def cb(sub):
        value=sum(m.y[(s,t),(i,j)]+m.z[(s,t),(i,j)] for (i,j) in m.Node[s].FS)
        return value==1
    
    m.Sub[(s,t)].sub.cb=Constraint(rule=cb)
    
    def cc(sub,i):
        value=sum(m.y[(s,t),(i,j)]+m.z[(s,t),(i,j)] for (i,j) in m.Node[i].FS)
        value=value - sum(m.p[(j,i)]*m.y[(s,t),(j,i)]+m.q[(j,i)]*m.z[(s,t),(j,i)] for (j,i) in m.Node[i].RS)
        return value == 0
    
    m.Sub[(s,t)].sub.cc=Constraint(m.nodeset-m.Sub[(s,t)].sub.stset,rule=cc)
    
    def cd(sub):
        value = sum (m.p[(j,i)]*m.y[(s,t),(j,i)]+m.q[(j,i)]*m.z[(s,t),(j,i)] for (j,i) in m.Node[t].RS)
        return m.h[(s,t)]-value == 0
    m.Sub[(s,t)].sub.cd=Constraint(rule=cd)

    def ce(sub,i,j):
        return m.y[(s,t),(i,j)] <= 1- m.x[(i,j)]
    
    m.Sub[(s,t)].sub.ce=Constraint(m.Edgeset,rule=ce)
    
    def cf(sub,i,j):
        return m.z[(s,t),(i,j)] <= m.x[(i,j)]
    
    m.Sub[(s,t)].sub.cf=Constraint(m.Edgeset,rule=cf)

    
m.pprint()

kwargs={'use_dual_objective':False}
solver = SolverFactory('pao.bilevel.ld')  
solver.options.solver = 'gurobi'  
results = solver.solve(m, **kwargs, tee=True)  
m.y.pprint()
m.z.pprint()
'''
for (i,j) in m.Edgeset:
    print(f'x_{(i,j)}={m.x[(i,j)].value}')
    for (s,t) in m.OmegaSet:
        print(f'y_{(i,j)}={m.Sub[(s,t)].sub.y[(i,j)].value}, z_{(i,j)}={m.Sub[(s,t)].sub.z[(i,j)].value}')
'''    
