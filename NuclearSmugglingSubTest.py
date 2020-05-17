# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 19:23:17 2020

@author: sheifazera
"""

import time

from pyomo.environ import *
from pao.bilevel import *

from NuclearSmugglingParameters20e import E, FS, RS, n, e, b, c, p, q, P, ns, nt, W

s=5;
t=8;


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
    
x={}
x[(3,7)]=0
x[(4,6)]=0
x[(5,4)]=0
x[(5,8)]=1
x[(9,12)]=0
x[(11,9)]=0
x[(2,1)]=0
x[(3,1)]=0
x[(3,2)]=0
x[(3,6)]=0
x[(4,3)]=0
x[(4,5)]=0
x[(5,2)]=0
x[(5,6)]=0
x[(6,4)]=0
x[(6,7)]=0
x[(6,8)]=0
x[(6,9)]=0
x[(9,11)]=0
x[(10,6)]=0
x[(11,10)]=0
x[(11,12)]=1

m.p=Param(m.Edgeset,initialize=p)
m.q=Param(m.Edgeset,initialize=q)
m.x=Param(m.Edgeset,initialize=x)
m.h=Var(m.OmegaSet,within=NonNegativeReals)
m.y=Var(m.OmegaSet,m.Edgeset,within=NonNegativeReals)
m.z=Var(m.OmegaSet,m.Edgeset,within=NonNegativeReals)



m.stset=Set(initialize=[s,t])
def obj(sub):
    value=m.h[(s,t)]
    return value
m.obj=Objective(rule=obj, sense=maximize)
    
def cb(sub):
    value=sum(m.y[(s,t),(i,j)]+m.z[(s,t),(i,j)] for (i,j) in m.Node[s].FS)
    return value==1
    
m.cb=Constraint(rule=cb)
    
def cc(sub,i):
    value=sum(m.y[(s,t),(i,j)]+m.z[(s,t),(i,j)] for (i,j) in m.Node[i].FS)
    value=value - sum(m.p[(j,i)]*m.y[(s,t),(j,i)]+m.q[(j,i)]*m.z[(s,t),(j,i)] for (j,i) in m.Node[i].RS)
    return value == 0
    
m.cc=Constraint(m.nodeset-m.stset,rule=cc)
    
def cd(sub):
    value = sum (m.p[(j,i)]*m.y[(s,t),(j,i)]+m.q[(j,i)]*m.z[(s,t),(j,i)] for (j,i) in m.Node[t].RS)
    return m.h[(s,t)]-value == 0
m.cd=Constraint(rule=cd)

def ce(sub,i,j):
    return m.y[(s,t),(i,j)] <= 1- m.x[(i,j)]
    
m.ce=Constraint(m.Edgeset,rule=ce)
    
def cf(sub,i,j):
    return m.z[(s,t),(i,j)] <= m.x[(i,j)]
    
m.cf=Constraint(m.Edgeset,rule=cf)


solver=SolverFactory('gurobi')
results=solver.solve(m,tee=True)

for (i,j) in m.Edgeset:
    print(f'y_{(i,j)}={m.y[(s,t),(i,j)].value}, z_{(i,j)}={m.z[(s,t),(i,j)].value}')
