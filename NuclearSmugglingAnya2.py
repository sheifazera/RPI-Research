# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 12:34:08 2020

@author: sheifazera
"""

import time

from pyomo.environ import *
from pao.bilevel import *

from NuclearSmugglingParametersAnya2 import E, FS, RS, n, e, b, c, p, q, P, ns, nt, W

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

def Obj(m):
    value=sum(m.P[(i,j)]*m.h[(i,j)] for (i,j) in m.OmegaSet) #h[i,j] is the objective of subproblem (i,j)
    return value
m.Obj=Objective(rule=Obj,sense=minimize)

def C(m):
    value=sum(m.c[(i,j)]*m.x[(i,j)] for (i,j) in m.Edgeset)
    return value <= b
m.C=Constraint(rule=C)


m.Sub=Block(m.OmegaSet)
for (s,t) in m.OmegaSet:
    m.Sub[(s,t)].sub=SubModel(fixed=m.x)
    m.Sub[(s,t)].sub.stset=Set(initialize=[s,t])
    m.Sub[(s,t)].sub.y=Var(m.Edgeset,within=NonNegativeReals, bounds=(0,1))
    m.Sub[(s,t)].sub.z=Var(m.Edgeset,within=NonNegativeReals, bounds=(0,1))
    def obj(sub):
        value=m.h[(s,t)]
        return value
    m.Sub[(s,t)].sub.obj=Objective(rule=obj, sense=maximize)

    def cb(sub):
        value=sum(m.Sub[(s,t)].sub.y[(i,j)]+m.Sub[(s,t)].sub.z[(i,j)] for (i,j) in m.Node[s].FS)
        return value==1
    m.Sub[(s,t)].sub.cb=Constraint(rule=cb)
    
    def cc(sub,i):
        value=sum(m.Sub[(s,t)].sub.y[i,j]+m.Sub[(s,t)].sub.z[(i,j)] for (i,j) in m.Node[i].FS)
        value=value - sum(m.p[(j,i)]*m.Sub[(s,t)].sub.y[(j,i)]+m.q[(j,i)]*m.Sub[(s,t)].sub.z[(j,i)] for (j,i) in m.Node[i].RS)
        return value == 0
    m.Sub[(s,t)].sub.cc=Constraint(m.nodeset-m.Sub[(s,t)].sub.stset,rule=cc)
    
    def cd(sub):
        value = sum (m.p[(j,i)]*m.Sub[(s,t)].sub.y[(j,i)]+m.q[(j,i)]*m.Sub[(s,t)].sub.z[(j,i)] for (j,i) in m.Node[t].RS)
        return m.h[(s,t)]-value == 0
    m.Sub[(s,t)].sub.cd=Constraint(rule=cd)
    
    def ce(sub,i,j):
        return m.Sub[(s,t)].sub.y[(i,j)] <= 1- m.x[(i,j)]
    m.Sub[(s,t)].sub.ce=Constraint(m.Edgeset,rule=ce)
    
    def cf(sub,i,j):
        return m.Sub[(s,t)].sub.z[(i,j)] <= m.x[(i,j)]
    m.Sub[(s,t)].sub.cf=Constraint(m.Edgeset,rule=cf)

#m.pprint()
weights = {m.Sub[k].name+'.sub':v for k,v in m.P.items()}
kwargs={'subproblem_objective_weights':weights, 'use_dual_objective':False}
solver = SolverFactory('pao.bilevel.stochastic_ld')
solver.options.solver = 'gurobi'  
results = solver.solve(m, **kwargs, tee=True)
for (i,j) in m.Edgeset:
    print(f'x_{(i,j)}={m.x[(i,j)].value}')
for (s,t) in m.OmegaSet:
    print('Subproblem ({},{})'.format(s,t))
    for (i, j) in m.Edgeset:
        print(f'y_{(i,j)}={m.Sub[(s,t)].sub.y[(i,j)].value}, z_{(i,j)}={m.Sub[(s,t)].sub.z[(i,j)].value}')
print('done')
