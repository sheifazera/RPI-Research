# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 11:41:37 2020

@author: sheifazera
"""

import numpy as np
from pyomo.environ import *
'''
    min cR*xu + cZ*yu + dR*xl0 + dZ*yl0 
    s.t. AR*xu + AZ*yu + BR*xl0 + BZ* yl0 <= r
'''

mU = 2
nL= 2 

nR = 1

nZ = 1

mR = 1

mZ = 1
 
AR_array = np.array([0, 6]).reshape(2,1).tolist()
AR={}
for i in range(0,2):
    for j in range(0,1):
        AR[(i+1,j+1)]=AR_array[i][j]
        
        
AZ_array = np.array([7,9]).reshape(2,1).tolist()
AZ={}
for i in range(0,2):
    for j in range(0,1):
        AZ[(i+1,j+1)]=AZ_array[i][j]

BR_array = np.array([5, 10]).reshape(2,1).tolist()
BR={}
for i in range(0,2):
    for j in range(0,1):
        BR[(i+1,j+1)]=BR_array[i][j]

BZ_array = np.array([7,2]).reshape(2,1).tolist()
BZ={}
for i in range(0,2):
    for j in range(0,1):
        BZ[(i+1,j+1)]=BZ_array[i][j]

r_array = np.array([62, 117]).reshape(2,1).tolist()
r={}
for i in range(0,2):
    r[i+1]=r_array[i][0]    
    
cR_array = np.array([20]).reshape(1,1).tolist()
cR={}
for i in range(0,1):
    cR[i+1]=cR_array[i][0]
 
cZ_array =  np.array([-38]).reshape(1,1).tolist()
cZ={}
for i in range(0,1):
    cZ[i+1]=cZ_array[i][0]
dR_array = np.array([1]).reshape(1,1).tolist()
dR={}
for i in range(0,1):
    dR[i+1]=dR_array[i][0]
    
dZ_array =  np.array([42]).reshape(1,1).tolist()
dZ={}
for i in range(0,1):
    dZ[i+1]=dZ_array[i][0]
    
    
#Turn it into a Pyomo model:


Parent=ConcreteModel()

Parent.mUset=RangeSet(1,mU)
Parent.nLset=RangeSet(1,nL)
Parent.nRset=RangeSet(1,nR)
Parent.nZset=RangeSet(1,nZ)
Parent.mRset=RangeSet(1,mR)
Parent.mZset=RangeSet(1,mZ)

Parent.r=Param(Parent.mUset,initialize=r,default=0,mutable=True)
Parent.cR=Param(Parent.mRset,initialize=cR,default=0,mutable=True)
Parent.cZ=Param(Parent.mZset,initialize=cZ,default=0,mutable=True)
Parent.dR=Param(Parent.nRset,initialize=dR,default=0,mutable=True)
Parent.dZ=Param(Parent.nZset,initialize=dZ,default=0,mutable=True)


Parent.AR=Param(Parent.mUset,Parent.mRset,initialize=AR, default=0,mutable=True)
Parent.AZ=Param(Parent.mUset,Parent.mZset,initialize=AZ, default=0,mutable=True)
Parent.BR=Param(Parent.mUset,Parent.nRset,initialize=BR, default=0,mutable=True)
Parent.BZ=Param(Parent.mUset,Parent.nZset,initialize=BZ, default=0,mutable=True)

Parent.zero=Param(initialize=0, mutable=True) 


Parent.Master=Block()

Parent.Master.Y=Param(Any,within=NonNegativeIntegers,mutable=True) #Growing Parameter
Parent.Master.xu=Var(Parent.mRset,within=NonNegativeReals)
Parent.Master.yu=Var(Parent.mZset,within=NonNegativeIntegers)
Parent.Master.xl0=Var(Parent.nRset,within=NonNegativeReals)
Parent.Master.yl0=Var(Parent.nZset,within=NonNegativeIntegers)


def Master_obj(Master):
    value=(sum(Parent.cR[j]*Parent.Master.xu[j] for j in Parent.mRset)+
           sum(Parent.cZ[j]*Parent.Master.yu[j] for j in Parent.mZset)+
           sum(Parent.dR[j]*Parent.Master.xl0[j] for j in Parent.nRset)+
           sum(Parent.dZ[j]*Parent.Master.yl0[j] for j in Parent.nZset))
    return value

Parent.Master.Theta_star=Objective(rule=Master_obj,sense=minimize)
    
def Master_c1(Master,i):
    value=(sum(Parent.AR[(i,j)]*Parent.Master.xu[j] for j in Parent.mRset)+
           sum(Parent.AZ[(i,j)]*Parent.Master.yu[j] for j in Parent.mZset)+
           sum(Parent.BR[(i,j)]*Parent.Master.xl0[j] for j in Parent.nRset)+
           sum(Parent.BZ[(i,j)]*Parent.Master.yl0[j] for j in Parent.nZset))
    return value - Parent.r[i] <= Parent.zero
Parent.Master.c1=Constraint(Parent.mUset,rule=Master_c1) #(12)    