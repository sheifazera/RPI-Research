# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 16:21:06 2020

@author: sheifazera
"""

from pyomo.environ import *
from pyomo.gdp import *

opt=SolverFactory("gurobi")
bigm=TransformationFactory('gdp.bigm')



m=ConcreteModel()
m.x=Var(within=NonNegativeReals)
m.y=Var(within=NonNegativeReals)
m.z=Var(within=NonNegativeReals)

m.a=Param(initialize=2,mutable=True)
m.b=Param(initialize=1,mutable=True)
m.c=Param(initialize=1,mutable=True)

m.d=Param(initialize=3, mutable=True)

# Primal Problem
def obj(m):
    value=m.a*m.x+m.b*m.y+m.c*m.z
    return value

m.obj=Objective(rule=obj,sense=minimize)

def c1(m):
    value=m.y*m.y+m.z*m.z
    return value <= m.x*m.x

m.c1=Constraint(rule=c1)

def c2(m):
    value=m.x+m.y+m.z
    return value==m.d

m.c2=Constraint(rule=c2)


Presults=opt.solve(m)
print(f'x={m.x.value},y={m.y.value},z={m.z.value}')

#Dual Problem

d=ConcreteModel()
d.a=Param(initialize=2,mutable=True)
d.b=Param(initialize=1,mutable=True)
d.c=Param(initialize=1,mutable=True)

d.d=Param(initialize=3, mutable=True)
d.eta1=Var(within=NonNegativeReals)
d.eta2=Var(within=Reals)
d.eta3=Var(within=Reals)

d.pi=Var(within=Reals)
def objdual(d):
    value=d.d*d.pi
    return value

d.obj=Objective(rule=objdual,sense=maximize)

def c1dual(d):
    value=d.pi+d.eta1
    return value ==d.a

d.c1=Constraint(rule=c1dual)

def c2dual(d):
    value=d.pi+d.eta2
    return value ==d.b

d.c2=Constraint(rule=c2dual)

def c3dual(d):
    value=d.pi+d.eta3
    return value ==d.c

d.c3=Constraint(rule=c3dual)

def c4dual(d):
    value=d.eta2*d.eta2+d.eta3*d.eta3
    return value <= d.eta1*d.eta1

d.c4=Constraint(rule=c4dual)

#d.pprint()
Dresults=opt.solve(d)
print(f'xpi={d.pi.value},eta1={d.eta1.value},eta2={d.eta2.value},eta3={d.eta3.value}')

#Dual Optimality 
o=ConcreteModel()

o.x=Var(within=NonNegativeReals)
o.y=Var(within=NonNegativeReals)
o.z=Var(within=NonNegativeReals)

o.a=Param(initialize=2,mutable=True)
o.b=Param(initialize=1,mutable=True)
o.c=Param(initialize=1,mutable=True)

o.d=Param(initialize=3, mutable=True)
o.eta1=Var(within=NonNegativeReals)
o.eta2=Var(within=Reals)
o.eta3=Var(within=Reals)
o.pi=Var(within=Reals)

o.c1d=Constraint(rule=c1dual)

o.c2d=Constraint(rule=c2dual)

o.c3d=Constraint(rule=c3dual)

o.c4d=Constraint(rule=c4dual)

o.c1p=Constraint(rule=c1)

o.c2p=Constraint(rule=c2)

def comp1(o):
    value=o.eta1*o.y+o.x*o.eta2
    return value==0
o.comp1=Constraint(rule=comp1)

def comp2(o):
    value=o.eta1*o.z+o.x*o.eta3
    return value==0
o.comp2=Constraint(rule=comp2)

def OptObj(o):
    value=0
    return value

o.obj=Objective(rule=OptObj,sense=maximize)
bigm.apply_to(o)
o.pprint()


Optresults=opt.solve(o)
print(f'x={o.x.value},y={o.y.value},z={o.z.value}')
print(f'xpi={o.pi.value},eta1={o.eta1.value},eta2={o.eta2.value},eta3={o.eta3.value}')

