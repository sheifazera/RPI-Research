# -*- coding: utf-8 -*-
'''
Created on Thu Oct 17 13:55:54 2019

@author: She'ifa
'''

from pyomo.environ import *
from pyomo.gdp import *
from pyomo.mpec import *

infinity = float('inf')
opt=SolverFactory("gurobi")
bigm = TransformationFactory('gdp.bigm')

from TestExample import c,A, b, m, n

m=ConcreteModel()
m.mSet=RangeSet(1,m)
m.nSet=RangeSet(1,n)

m.c=Param(m.nSet,initialize=c)
m.A=Param(m.mSet,m.nSet,initialize=A, mutable=true)
m.b=Param(m.mSet,initialize=b,mutable=true)
m.x=Var(m.nSet,within=Reals)
m.y=Var(within=Reals)

def obj(m):
    value=sum(m.c[i]*m.x[i] for i in m.nSet)
    return value
m.Obj=Objective(rule=obj,sense=minimize) #min c^Tx


def c1(m):
    value=0
    for i in m.nSet:
        value=value+m.x[i]*m.x[i]
    return value<=m.y*m.y
m.c1=Constraint(rule=c1)

def c2(m):
    value=sum(m.x[i] for i in m.nSet) +m.y
    return value <= 1
m.c2=Constraint(rule=c2)

results=opt.solve(m) 

    