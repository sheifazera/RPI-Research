# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 13:12:22 2020

@author: She'ifa
"""

from pyomo.environ import *
from pyomo.gdp import *
from pyomo.mpec import *

infinity = float('inf')
opt=SolverFactory("gurobi")
bigm = TransformationFactory('gdp.bigm')

from ToyExample2Uncertain import mU, nL, nR, nZ, mR, mZ,  PR, PZ, s, QR, QZ, wR, wZ, r, NR, NZ

Parent=ConcreteModel()

Parent.mUset=RangeSet(1,mU)
Parent.nLset=RangeSet(1,nL)
Parent.nRset=RangeSet(1,nR)
Parent.nZset=RangeSet(1,nZ)
Parent.mRset=RangeSet(1,mR)
Parent.mZset=RangeSet(1,mZ)


Parent.alpha=Var(Parent.nLset,within=Reals)
Parent.beta=Var(Parent.nRset,within=Reals)
Parent.gamma=Var(Parent.nZset, within=Reals)
Parent.tau=Var(within=NonNegativeReals)

Parent.s=Param(Parent.nLset,initialize=s,default=0,mutable=True)
Parent.PR=Param(Parent.nLset,Parent.nRset,initialize=PR, default=0,mutable=True)
Parent.PZ=Param(Parent.nLset,Parent.nZset,initialize=PZ, default=0,mutable=True)
Parent.QR=Param(Parent.nLset,Parent.mRset,initialize=QR, default=0,mutable=True)
Parent.QZ=Param(Parent.nLset,Parent.mZset,initialize=QZ, default=0,mutable=True)
Parent.NR=Param(Parent.nRset,Parent.nRset,initialize=NR,mutable=True)
Parent.NZ=Param(Parent.nZset,Parent.nZset,initialize=NZ,mutable=True)
Parent.wR=Param(Parent.nRset,initialize=wR,default=0,mutable=True)

Parent.xu=Param(Parent.mRset,initialize=1,mutable=True)
Parent.yu=Param(Parent.mZset,initialize=1,mutable=True)
Parent.yl=Param(Parent.nZset,initialize=1,mutable=True)
Parent.w=Var(within=NonNegativeReals)
Parent.nx=Var(Parent.nRset,within=Reals)
Parent.ny=Var(Parent.nZset,within=Reals)
Parent.NZY=Param(Parent.nZset, initialize=1, mutable=True)
for i in Parent.nZset:
    Parent.NZY[i]=sum(Parent.NZ[i,j]*Parent.yl[j] for j in Parent.nZset)


def Obj(Parent):
    value= sum((Parent.s[i]-sum(Parent.QR[i,j]*Parent.xu[j] for j in Parent.mRset)-
                sum(Parent.QZ[i,j]*Parent.yu[j] for j in Parent.mZset)-
                sum(Parent.PZ[i,j]*Parent.yl[j] for j in Parent.nZset))*Parent.alpha[i] for i in Parent.nLset)

    value=value-sum(Parent.NZY[i]*Parent.gamma[i] for i in Parent.nZset)
    return value

Parent.Obj=Objective(rule=Obj,sense=minimize)
'''
def c1(Parent):
    value=
    return value
Parent.c1=Constraint(rule=c1)
'''
def c2(Parent):
    value=sum(Parent.alpha[i] for i in Parent.nLset)-Parent.tau
    return value==0
Parent.c2=Constraint(rule=c2)

def c3(Parent):
    value=sum(Parent.beta[i]*Parent.beta[i] for i in Parent.nRset)
    value=value+sum(Parent.gamma[i]*Parent.gamma[i] for i in Parent.nZset)
    return value <= Parent.tau*Parent.tau
Parent.c3=Constraint(rule=c3)

def c4a(Parent,i):
    value=Parent.tau*Parent.ny[i] - Parent.w*Parent.gamma[i]
    return value == 0
Parent.c4a=Constraint(Parent.nZset,rule=c4a)

def c4b(Parent,i):
    value= Parent.tau*Parent.nx[i] - Parent.w*Parent.beta[i]
    return value == 0
Parent.c4b=Constraint(Parent.nRset,rule=c4b)

Parent.pprint()


results=opt.solve(Parent)
