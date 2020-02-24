# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 10:33:00 2020

@author: She'ifa
NORVEP: Besancon, et al. 2019
"""

from pyomo.environ import *
from pyomo.gdp import *
from pyomo.mpec import *
#from pao.bilevel import *
#from pao.duality import * 
import cdd

infinity = float('inf')
opt=SolverFactory("gurobi")
bigm = TransformationFactory('gdp.bigm')
from NORVEPTestExample import cx, cy, G, H, q,  d, A, B, b, delta, nl, ml,nu,mu
#The upper level variable is x in $R_+^{n_u}
#The lower level variable is v in R_+^{n_l}
#cx is nu
#cy is nl
#G is mu x nu
#H is mu x nl
#q is mu
#d is nl
#A is ml x nu
#B is ml x nl
#b is ml

m=ConcreteModel()
m.nlset=RangeSet(1,nl)
m.nuset=RangeSet(1,nu)
m.mlset=RangeSet(1,ml)
m.muset=RangeSet(1,mu)
m.x=Var(m.nuset,within=NonNegativeReals)
m.v=Var(m.nlset,within=NonNegativeReals)
m.cx=Param(m.nuset,initialize=cx,mutable=True)
m.cy=Param(m.nlset,initialize=cy,mutable=True)
m.A=Param(m.mlset,m.nuset,initialize=A,mutable=True)
m.B=Param(m.mlset,m.nlset,initialize=B,mutable=True)
m.b=Param(m.mlset,initialize=b,mutable=True)
m.G=Param(m.muset,m.nuset,initialize=G,mutable=True)
m.H=Param(m.muset,m.nlset,initialize=H,mutable=True)
m.q=Param(m.muset,initialize=q,mutable=True)
m.d=Param(m.nlset,initialize=d,mutable=True)


def Obj(m):
    value=sum(m.cx[i]*m.x[i] for i in m.nuset)
    vaule=value + sum(m.cy[i]*m.v[i] for i in m.nlset)
    return value
m.Obj=Objective(rule=Obj,sense=minimize)
def c1(m,i):
    value=sum(m.G[i,j]*m.x[j] for j in m.nuset)
    value=value + sum(m.H[i,j]*m.v[j] for j in m.nlset)
    return value <= m.q[i]
m.c1=Constraint(m.muset,rule=c1)

#Lower Level
m.LL=SubModel(fixed=m.x)
def LObj(m):
    value=sum(m.d[i]*m.v[i] for i in m.nlset)
    return value
m.LL.LObj=Objective(rule=LObj,sense=maximize)
def c2(m,i):
    value=sum(m.A[i,j]*m.x[j] for j in m.nuset)
    value=value + sum(m.B[i,j]*m.v[j] for j in m.nlset)
    return value <= m.b[i]
m.LL.c2=Constraint(m.mlset,rule=c2)

#Highpoint Problem
'''
highpoint = TransformationFactory('pao.bilevel.highpoint')
HP=highpoint.apply_to(m)
resultsHP=opt.solve(HP)
if resultsHP.solver.termination_condition==TerminationCondition.infeasible or resultsHP.solver.termination_condition==TerminationCondition.infeasibleOrUnbounded: 
    raise RuntimeError("Highpoint Problem infeasible")
'''  
#Highpoint Problem (Code Option 2)
highpoint = TransformationFactory('pao.bilevel.highpoint')
highpoint.apply_to(m)

solver = SolverFactory('gurobi')
for c in instance.component_objects(Block, descend_into=False):
    if '_hp' in c.name:
        c.activate()
        results = solver.solve(c, tee=True, keepfiles=True)
        c.deactivate()
 
if results.solver.termination_condition==TerminationCondition.infeasible or results.solver.termination_condition==TerminationCondition.infeasibleOrUnbounded: 
    raise RuntimeError("Highpoint Problem infeasible")

#Optimistic Bilevel 
''' 
opt2 = SolverFactory(pao.bilevel.blp_local)
opt2.options.solver = 'gurobi'
resultsOB=opt2.solve(m)

if resultsOB.solver.termination_condition==TerminationCondition.infeasible or resultsOB.solver.termination_condition==TerminationCondition.infeasibleOrUnbounded: 
    raise RuntimeError("Optimistic Bilevel Problem infeasible")
'''

#Adversarial
'''
m.alpha=Var(m.mlset,within=NonNegativeReals)
m.beta=Var(within=NonNegativeReals)
m.Vertices=Parameter(Any)
m.Adversarial=Block(m.muset)
def c3(m,i,k):
    value=sum(m.B[j,i]*m.alpha[j] for j in m.mlset) + m.b[i]*m.beta
    value<=m.H[i,k]
    return value

#Get Bd=[B^T|d] matrix for all adversarial problems using numpy

for k in m.muset:
    #Get Hk as numpy array
    Hk=
    #Concatenate [-Hk,-[B^T|d]]
    mat=cdd.Matrix(hstack((Hk,Bd)),number_type='fraction')
    mat.rep_type=cdd.RepType.INEQUALITY
    poly=cdd.Polyhedron(mat)
    ext=poly.get_generators()
    #Get the extreme points
    m.Adversarial[k].alpha=Variable(m.mlset,within=NonNegativeReals)
    m.Adversarial[k].beta=Variable(within=NonNegativeReals) #Scalar
    m.Adversarial[k].c3=Constraint(m.muset,rule=c3)
    results=opt.solve(m.Adversarial[k])
    for l in #vertices of the k-th polytope
    m.Vertices[k,l]=results.alpha.value #Vertex l of the k-th polytope
'''    
    
    

#Extended Aggragated Near-Optimal Problem

