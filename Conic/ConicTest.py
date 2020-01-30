# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 12:26:35 2019

@author: She'ifa
"""

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

from TestExample import c,A, b, m, n, P

model=ConcreteModel()
model.mSet=RangeSet(1,m)
model.nSet=RangeSet(1,n)

model.c=Param(model.nSet,initialize=c)
model.A=Param(model.mSet,model.nSet,initialize=A, mutable=True)
model.b=Param(model.mSet,initialize=b,mutable=True)
model.x=Var(model.nSet,within=NonNegativeReals)
model.y=Var(within=NonNegativeReals)
model.P=Param(model.mSet,model.nSet,initialize=P,mutable=True)
model.Px=Var(model.mSet,within=NonNegativeReals) #If P is not positive, then this is not necessarily positive
model.t=Var(within=NonNegativeReals)

def obj(model):
    value=sum(model.c[i]*model.x[i] for i in model.nSet)
    return value
model.Obj=Objective(rule=obj,sense=minimize) #min c^Tx

def c1(model,j):
    value=sum(model.P[j,i]*model.x[i] for i in model.nSet)
    return value==model.Px[j]
model.c1=Constraint(model.mSet,rule=c1) #Creating a Variable Px=P^Tx

def c2(model):
    value=sum(model.Px[i]*model.Px[i] for i in model.mSet)
    return value <= model.t*model.t
model.c2=Constraint(model.mSet,rule=c2)

def c3(model,j):
    value=sum(model.A[j,i]*model.x[i] for i in model.nSet)+model.t
    return value <= model.b[j]
model.c3=Constraint(model.mSet,rule=c3)

results=opt.solve(model) 

    