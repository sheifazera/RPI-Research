# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 11:59:00 2020

@author: sheifazera
"""

from pyomo.environ import * 

opt=SolverFactory("gurobi")

m=ConcreteModel()
m.x=Var(within=NonNegativeReals)
m.a=Param(initialize=1)

def c1(m):
    value=m.a*m.x
    return value <= 1
m.c1=Constraint(rule=c1)

def obj(m):
    value=m.x
    return value
m.obj=Objective(rule=obj,sense=minimize)

opt.solve(m)
print(f'x={m.x.value}')