# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 20:18:30 2020

@author: sheifazera
"""

'''
Implementation of the algorithm presented in "A projection-based reformulation
and decomposition algorithm for global optimization of a class of mixed integer 
bilevel linear programs" by Dajun Yue, Jiyao Gao, Bo Zeng, Fengqi You

Implemented May-August 2019 by She'ifa Punla-Green at Sandia National Labs

This algorithm seeks to solve the following bilevel MILP:
    min cR*xu + cZ*yu + dR*xl0 + dZ*yl0 
    s.t. AR*xu + AZ*yu + BR*xl0 + BZ* yl0 <= r
     (xl0,yl0) in argmax {wR*xl+wZ*yl: PR*xl+PZ*yl<=s-QR*xu-QZ*yu}
'''
import time
import sys 

from pyomo.environ import *
from pyomo.gdp import *
from pyomo.mpec import *


solvername='SCIPAMPL'

solverpath_folder='C:\\Users\spunlag\Documents\SCIP' 

solverpath_exe='C:\\Users\spunlag\Documents\SCIP\scipampl'

sys.path.append(solverpath_folder)

SCIP=SolverFactory(solvername,executable=solverpath_exe)
#SCIP=SolverFactory('SCIPAMPL')


infinity = float('inf')
IPOPT=SolverFactory("ipopt")
GUROBI=SolverFactory("gurobi")
bigm = TransformationFactory('gdp.bigm')

'''
Parameter Import

Notation:
m is for upper level
n is for lower level
x is continuous (R)
y is discrete (Z)

coefficient vectors and matrices should be dictionaries (with tuples for matrices)
variable size should be floats
'''
from ToyExample2Uncertain import mU, nL, nR, nZ, mR, mZ, AR, AZ, BR, BZ, cR, cZ, dR, dZ, PR, PZ, s, QR, QZ, wR, wZ, r, zx, zy, NZ, NR
'''
mU number of upper level constraints
nL number of lower level constraints
nR number of lower level continuous variables
nZ number of lower level integer variables
mR number of upper level continuous variables
mZ number of upper level integer variables
AR constraint matrix for upper level problem, upper level continuous variables
AZ constraint matrix for upper level problem, upper level integer variables
BR constraint matrix for upper level problem, lower level continuous variables
BZ constraint matrix for upper level problem, lower level integer variables
r  RHS vector for upper level constraint
cR coefficient vector for upper level objective, upper level continuous variables
cZ coefficient vector for upper level objective, upper level integer variables
dR coefficient vector for upper level objective, lower level continuous variables
dZ coefficient vector for upper level objective, lower level integer variables
PR constraint matrix for lower level problem, lower level continuous variables
PZ constraint matrix for lower level problem, lower level integer variables
QR constraint matrix for lower level problem, upper level continuous variables
QZ constraint matrix for lower level problem, upper level integer variables
s  RHS vector for lower level constraint
wR coefficient vector for the lower level objective, lower level continuous variables
wZ coefficient vector for the lower level objective, lower level integer variables
'''

t = time.time()

#These parameters can be changed for your specific problem
epsilon= 1e-4 #For use in disjunction approximation
xi=0 #tolerance for UB-LB to claim convergence
maxit=5 #Maximum number of iterations
M=1e6 #upper bound on variables

LB=-infinity
UB=infinity
k=0
 
flag=0

#Master problem, subproblem 1, and subproblem 2 are all blocks on a Parent concrete model
#so that parameters that are present in all three problems can be shared
#They are mutable parameters on a parent model rather than just python dictionaries so that 
#Pyomo constraints that would be trivially true such as 0<=0 do not cause an error

Parent=ConcreteModel()

Parent.mUset=RangeSet(1,mU)
Parent.nLset=RangeSet(1,nL)
Parent.nRset=RangeSet(1,nR)
Parent.nZset=RangeSet(1,nZ)
Parent.mRset=RangeSet(1,mR)
Parent.mZset=RangeSet(1,mZ)
Parent.zxset=RangeSet(1,zx)
Parent.zyset=RangeSet(1,zy)

Parent.s=Param(Parent.nLset,initialize=s,default=0,mutable=True)
Parent.r=Param(Parent.mUset,initialize=r,default=0,mutable=True)
Parent.cR=Param(Parent.mRset,initialize=cR,default=0,mutable=True)
Parent.cZ=Param(Parent.mZset,initialize=cZ,default=0,mutable=True)
Parent.dR=Param(Parent.nRset,initialize=dR,default=0,mutable=True)
Parent.dZ=Param(Parent.nZset,initialize=dZ,default=0,mutable=True)
Parent.wR=Param(Parent.nRset,initialize=wR,default=0,mutable=True)
Parent.wZ=Param(Parent.nZset,initialize=wZ,default=0,mutable=True)


Parent.AR=Param(Parent.mUset,Parent.mRset,initialize=AR, default=0,mutable=True)
Parent.AZ=Param(Parent.mUset,Parent.mZset,initialize=AZ, default=0,mutable=True)
Parent.BR=Param(Parent.mUset,Parent.nRset,initialize=BR, default=0,mutable=True)
Parent.BZ=Param(Parent.mUset,Parent.nZset,initialize=BZ, default=0,mutable=True)
Parent.PR=Param(Parent.nLset,Parent.nRset,initialize=PR, default=0,mutable=True)
Parent.PZ=Param(Parent.nLset,Parent.nZset,initialize=PZ, default=0,mutable=True)
Parent.QR=Param(Parent.nLset,Parent.mRset,initialize=QR, default=0,mutable=True)
Parent.QZ=Param(Parent.nLset,Parent.mZset,initialize=QZ, default=0,mutable=True)
Parent.NR=Param(Parent.zxset,Parent.nRset,initialize=NR, default=0,mutable=True)
Parent.NZ=Param(Parent.zyset,Parent.nZset,initialize=NZ, default=0,mutable=True)

Parent.zero=Param(initialize=0, mutable=True) 

'''
Master Problem: (P9)
    min cR*xu + cZ*yu + dR*xl0 + dZ*yl0 
    s.t. AR*xu + AZ*yu + BR*xl0 + BZ* yl0 <= r (12)
         PR*xl+PZ*yl<=s-QR*xu-QZ*yu (13)
         wR*xl0>= wR*xltilde (74)
         PR*xltilde <= s- QR*xu - QZ*yu-PZ*yl0 (75a)
         PR*pitilde>= wR (75b)
         xltilde \perp PR*pitilde - wR (76a)
         pitilde \perp s- QR*xu- QZ*yu - PR*xltilde (76b)
         PR*xlj-tj+PZ*yl<=s-QR*xu-QZ*yu for all 1<=j<=k (79)
         PR*lambdaj>=0 for all 1<=j<=k (83a)
         xlj \perp PR*lambdaj for all 1<=j<=k (83b)
         e-lambdaj >=0 for all 1<=j<=k (84a)
         tj \perp e - lambdaj for all 1<=j<=k (84b)
         lambdaj \perp s-QR*xu-QZ*yu-PR*xlj-PZ*ylj +tj for all 1<=j<=k (85)
         [e*tj=0] ==> [Constraint Block] for all 1<=j<=k (82)
         Constraint Block:
             wR*xl0 + wZ*yl0 >= wR*xlj + wZ*ylj
             PR*xlj <= s - QR*xu - QZ*yu - PZ*ylj
             PR*pij>= wR
             xlj \perp PR*pij-wR
             pij \perp s - QR*xu - QZ*yu - PZ*ylj - PR*xlj
'''

Parent.Master=Block()

Parent.Master.Y=Param(Any,within=NonNegativeIntegers,mutable=True) #Growing Parameter
Parent.Master.xu=Var(Parent.mRset,within=NonNegativeReals,bounds=(0,M))
Parent.Master.yu=Var(Parent.mZset,within=NonNegativeIntegers,bounds=(0,M))
Parent.Master.w=Var(within=NonNegativeReals,bounds=(0,M))
Parent.Master.xl0=Var(Parent.nRset,within=NonNegativeReals,bounds=(0,M))
Parent.Master.yl0=Var(Parent.nZset,within=NonNegativeIntegers,bounds=(0,M))
Parent.Master.nx=Var(Parent.zxset,within=Reals,bounds=(-M,M))
Parent.Master.ny=Var(Parent.zyset,within=Reals,bounds=(-M,M))

#Tilde Section
Parent.Master.xltilde=Var(Parent.nRset,within=NonNegativeReals,bounds=(0,M))
Parent.Master.wtilde=Var(within=NonNegativeReals,bounds=(0,M))
Parent.Master.nxtilde=Var(Parent.zxset,within=Reals,bounds=(-M,M))
Parent.Master.alphatilde=Var(Parent.nLset,within=NonNegativeReals,bounds=(0,M))
Parent.Master.betatilde=Var(Parent.zxset,within=Reals,bounds=(0,M))
Parent.Master.gammatilde=Var(Parent.zyset,within=Reals,bounds=(0,M))
Parent.Master.tautilde=Var(within=NonNegativeReals,bounds=(0,M))


#j section
Parent.Master.x=Var(Any,within=NonNegativeReals,dense=False,bounds=(0,M)) #Growing variable
Parent.Master.wj=Var(Any,within=NonNegativeReals,dense=False,bounds=(0,M))
Parent.Master.s=Var(Any,within=NonNegativeReals,dense=False,bounds=(0,M))
Parent.Master.nxj=Var(Any,within=Reals,dense=False,bounds=(-M,M))
Parent.Master.nyj=Var(Any,within=Reals,dense=False,bounds=(-M,M))
Parent.Master.alpha=Var(Any, within=NonNegativeReals,dense=False,bounds=(0,M)) #Growing variable
Parent.Master.beta=Var(Any,within=Reals,dense=False,bounds=(-M,M))
Parent.Master.gamma=Var(Any,within=Reals,dense=False,bounds=(-M,M))
Parent.Master.tau=Var(Any,within=NonNegativeReals,dense=False,bounds=(0,M))
Parent.Master.t=Var(Any,within=NonNegativeReals,dense=False,bounds=(0,M)) #Growing variable
Parent.Master.lam=Var(Any,within=NonNegativeReals,dense=False,bounds=(0,M)) #Growing variable
Parent.Master.mu=Var(Any,within=NonNegativeReals,dense=False,bounds=(0,M))
Parent.Master.eta=Var(Any,within=Reals,dense=False,bounds=(-M,M))
Parent.Master.nu=Var(Any,within=Reals,dense=False,bounds=(-M,M))

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

def Master_c2(Master,i): #ROBUST VERSION
    value=(sum(Parent.QR[(i,j)]*Parent.Master.xu[j] for j in Parent.mRset)+
           sum(Parent.QZ[(i,j)]*Parent.Master.yu[j] for j in Parent.mZset)+
           sum(Parent.PR[(i,j)]*Parent.Master.xl0[j] for j in Parent.nRset)+
           sum(Parent.PZ[(i,j)]*Parent.Master.yl0[j] for j in Parent.nZset))
    return value + Parent.Master.w - Parent.s[i] <= Parent.zero
Parent.Master.c2=Constraint(Parent.nLset,rule=Master_c2) #(13)

def Master_c2a(Master,i):
    alue=sum(Parent.NR[i,j]*Parent.Master.xl0[j] for j in Parent.nRset)
    return value==Parent.Master.nx[i]
Parent.Master.c2a=Constraint(Parent.zxset,rule=Master_c2a)

def Master_c2b(Master,i):
    value=sum(Parent.NZ[i,j]*Parent.Master.yl0[j] for j in Parent.nZset)
    return value==Parent.Master.ny[i]
Parent.Master.c2b=Constraint(Parent.zyset,rule=Master_c2b)

def Master_c2c(Master):
    value=sum(Parent.Master.nx[j]*Parent.Master.nx[j] for j in Parent.zxset )
    value=value+ sum(Parent.Master.ny[j]*Parent.Master.ny[j] for j in Parent.zyset)
    return value <= Parent.Master.w*Parent.Master.w
Parent.Master.c2c=Constraint(rule=Master_c2c)

def Master_c3(Master):
    l_value=sum(Parent.wR[i]*Parent.Master.xl0[i] for i in Parent.nRset)
    r_value=sum(Parent.wR[i]*Parent.Master.xltilde[i] for i in Parent.nRset)
    return l_value - r_value >= Parent.zero
Parent.Master.c3=Constraint(rule=Master_c3) #(74)

def Master_c4(Master,i): #ROBUST VERSION
    value=(sum(Parent.PR[(i,j)]*Parent.Master.xltilde[j] for j in Parent.nRset)+
           sum(Parent.PZ[(i,j)]*Parent.Master.yl0[j] for j in Parent.nZset)+
           sum(Parent.QZ[(i,j)]*Parent.Master.yu[j] for j in Parent.mZset)+
           sum(Parent.QR[(i,j)]*Parent.Master.xu[j] for j in Parent.mRset))
    return value + Parent.Master.wtilde - Parent.s[i] <= Parent.zero
Parent.Master.c4=Constraint(Parent.nLset,rule=Master_c4) #(75a)

def Master_c4a(Master,i):
    alue=sum(Parent.NR[i,j]*Parent.Master.xltilde[j] for j in Parent.nRset)
    return value==Parent.Master.nxtilde[i]
Parent.Master.c4a=Constraint(Parent.zxset,rule=Master_c4a)


def Master_c4c(Master):
    value=sum(Parent.Master.nxtilde[j]*Parent.Master.nxtilde[j] for j in Parent.zxset )
    value=value+ sum(Parent.Master.ny[j]*Parent.Master.ny[j] for j in Parent.zyset)
    return value <= Parent.Master.wtilde*Parent.Master.wtilde
Parent.Master.c4c=Constraint(rule=Master_c4c)

#Dual of lower level in terms of tilde variables

def Master_c5a(Master,j):
    value= sum(Parent.Master.alphatilde[i]*Parent.Master.PR[(j,i)] for i in Parent.nLset)
    value= value + sum(Parent.NR[j,i]*Parent.betatilde[i] for i in Parent.zxset)
    return value >= Parent.wR[j]
Parent.Master.c5a=Constraint(Parent.nRset,rule=Master_c5a)

def Master_c5b(Master):
    value=sum(Parent.Master.alphatilde[i] for i in Parent.nLset)
    return value - Parent.Master.tautilde ==0
Parent.Master.c5b=Constraint(rule=Master_c5b)

def Master_c5c(Master):
    value = sum( Parent.Master.betatilde[i]*Parent.Master.betatilde[i] for i in Parent.zxset)
    value= value + sum(Parent.Master.gammatilde[i]*Parent.Master.gammatilde[i] for i in Parent.zyset)
    return value <= Parent.Master.tautilde*Parent.Master.tautilde
Parent.Master.c5c=Constraint(rule=Master_c5c)

#Conic complementarity constraints
def Master_c5d(Master):
    value=Parent.Master.tautilde*Parent.Master.wtilde
    value=value - sum(Parent.Master.nxtilde[i]*Parent.Master.betatilde[i] for i in Parent.zxset)
    value=value - sum(Parent.Master.ny[i]*Parent.Master.gammatilde[i] for i in Parent.zyset)
    return value == Parent.zero
Parent.Master.c5d=Constraint(rule=Master_c5d)
    
def Master_c5e(Master,j):
    value=-Parent.Master.wtilde*Parent.Master.betatilde[j] - Parent.Master.tautilde*Parent.Master.nxtilde[j]
    return value == Parent.zero
Parent.Master.c5e=Constraint(Parent.zxset,rule=Master_c5e)
   
def Master_c5f(Master,j):
    value=Parent.Master.wtilde*Parent.Master.gammatilde[j] - Parent.Master.tautilde*Parent.Master.ny[j]
    return value == Parent.zero
Parent.Master.c5f=Constraint(Parent.zyset,rule=Master_c5f)

Parent.Master.CompBlock=Block()
Parent.Master.CompBlock.c6=ComplementarityList(rule=(complements(Parent.Master.xltilde[j] >= 0,
                                                       sum(Parent.Master.alphatilde[i]*Parent.Master.PR[(j,i)] for i in Parent.nLset) + sum(Parent.NR[j,i]*Parent.betatilde[i] for i in Parent.zxset)-
                                                       Parent.wR[j] <=0) for j in Parent.nRset))
#(76a)

#(76b)
Parent.Master.CompBlock.c7=ComplementarityList(rule=(complements(Parent.Master.alphatilde[j] >= 0,
                                                       (Parent.s[j]-sum(Parent.QR[(j,i)]*Parent.Master.xu[i] for i in Parent.mRset) -
                                                       sum(Parent.QZ[(j,i)]*Parent.Master.yu[i] for i in Parent.mZset)-
                                                       sum(Parent.PR[(j,i)]*Parent.Master.xltilde[i] for i in Parent.nRset)-
                                                       sum(Parent.PZ[(j,i)]*Parent.Master.yl0[i] for i in Parent.nZset)-Parent.Master.wtilde)>=0) for j in Parent.nLset))


Parent.Master.c_col = ConstraintList()
Parent.Master.CompBlock2=Block(Any)

Parent.Master.DisjunctionBlock=Block(Any)

def Master_add(Parent,k): #function for adding constraints on each iteration
    #(79)
    
    Parent.Master.CompBlock2[k].c_comp=ComplementarityList()
    for i in Parent.nLset:
        r_value= (sum(Parent.PR[(i,j)]*Parent.Master.x[(j,k)] for j in Parent.nRset)-
                  Parent.Master.t[(i,k)]+Parent.Master.wj[k])
        l_value= (Parent.s[i]-
                  sum(Parent.QR[(i,j)]*Parent.Master.xu[j] for j in Parent.mRset)-
                  sum(Parent.QZ[(i,j)]*Parent.Master.yu[j] for j in Parent.mZset)-
                  sum(Parent.PZ[(i,j)]*Parent.Master.Y[(j,k)] for j in Parent.nZset)) 
    
        Parent.Master.c_col.add(l_value-r_value >= Parent.zero)
    for i in Parent.zxset:
        r_value=sum(Parent.NR[(i,j)]*Parent.Master.x[(j,k)] for j in Parent.nRset)
        Parent.Master.c_col.add(r_value==Parent.Master.nxj[(i,k)])
    
    for i in Parent.zyset:
        r_value=sum(Parent.NZ[(i,j)]*Parent.Master.Y[(j,k)] for j in Parent.nZset)
        Parent.Master.c_col.add(r_value==Parent.Master.nyj[(i,k)])
    
    r_value=sum(Parent.Master.nxj[(j,k)]*Parent.Master.nxj[(j,k)] for j in Parent.zxset)
    r_value=r_value+sum(Parent.Master.nyj[(j,k)]*Parent.Master.nyj[(j,k)] for j in Parent.zyset)
    
    Parent.Master.c_col.add(r_value <= Parent.Master.wj[k]*Parent.Master.wj[k])
    
    #Optimality conditions for P9
    for i in Parent.nRset:  #Optimality for x
        Parent.Master.c_col.add(sum(Parent.Master.lam[i]*Parent.Master.PR[(j,i)] for i in Parent.Master.nLset) + sum(Parent.NR[j,i]*Parent.eta[i] for i in Parent.Master.zxset)<= Paremt.zero)
        Parent.Master.CompBlock2[k].c_comp.add(complements(Parent.Master.x[(i,k)]>=0, sum(Parent.Master.lam[i]*Parent.Master.PR[(j,i)] for i in Parent.Master.nLset) + sum(Parent.NR[j,i]*Parent.eta[i] for i in Parent.Master.zxset)>= 0))
                                                        #(83) 
     
    for i in Parent.nLset: #optimality for t #Already robust
        Parent.Master.c_col.add(1-Parent.Master.lam[(i,k)]>=Parent.zero) 
        Parent.Master.CompBlock2[k].c_comp.add(complements(1-Parent.Master.lam[(i,k)]>=0,Parent.Master.t[(i,k)]>=0)) #(84)
    
    for i in Parent.nLset: #Robust complementarity on first primal constraint
        Parent.Master.CompBlock2[k].c_comp.add(complements(Parent.Master.lam[(i,k)]>=0,(Parent.s[i]-
                                                     sum(Parent.QR[(i,j)]*Parent.Master.xu[j] for j in Parent.mRset)-
                                                     sum(Parent.QZ[(i,j)]*Parent.Master.yu[j] for j in Parent.mZset)-
                                                     sum(Parent.PZ[(i,j)]*Parent.Master.Y[(j,k)] for j in Parent.nZset)-
                                                     sum(Parent.PR[(i,j)]*Parent.Master.x[(j,k)] for j in Parent.nRset)+
                                                     Parent.Master.t[(i,k)]-Parent.Master.wj[k]>=0))) #(85)
        
    r_value=sum(Parent.Master.alpha[(i,k)] for j in Parent.nLset)-Parent.Master.mu[k]
    Parent.Master.c_col.add(r_value==0) #Dual on w
    
    r_value=sum(Parent.Master.eta[(j,k)]*Parent.Master.eta[(j,k)] for j in Parent.zxset)
    r_value=r_value+sum(Parent.Master.nu[(j,k)]*Parent.Master.nu[(j,k)] for j in Parent.zyset)
    Parent.Master.c_col.add(r_value <= Parent.Master.mu[k] * Parent.Master.mu[k]) # Dual Cone
    
    #Conic Complementarity
    
    r_value=Parent.Master.mu[k]*Parent.Master.wj[k]-sum(Parent.Master.eta[(j,k)]*Parent.Master.nxj[(j,k)] for j in Parent.zxset)-sum(Parent.Master.nu[(j,k)]*Parent.Master.nyj[(j,k)] for j in Parent.zyset)
    Parent.Master.c_col.add(r_value==0)
    
    for i in Parent.zxset:
        r_value=-Parent.Master.wj[k]*Parent.eta[(i,k)] + Parent.Master.mu[k]*Parent.nxj[(i,k)]
        Parent.Master.c_col.add(r_value==0)
    for i in Parent.zyset:     
        r_value=-Parent.Master.wj[k]*Parent.Master.nu[(i,k)] + Parent.Master.mu[k]*Parent.Master.nyj[(i,k)]
        Parent.Master.c_col.add(r_value==0)
        
    
    
    
    
    #(82) Disjunction
    Parent.Master.DisjunctionBlock[k].LH=Disjunct()
    Parent.Master.DisjunctionBlock[k].BLOCK=Disjunct()
    
    
    
    Parent.Master.DisjunctionBlock[k].LH.cons=Constraint(expr= sum(Parent.Master.t[(j,k)] for j in Parent.nLset) >= epsilon) 
    
    Parent.Master.DisjunctionBlock[k].BLOCK.cons=ConstraintList()
    
    l_value = (sum(Parent.wR[j]*Parent.Master.xl0[j] for j in Parent.nRset)+
               sum(Parent.wZ[j]*Parent.Master.yl0[j] for j in Parent.nZset))
    r_value = (sum(Parent.wR[j]*Parent.Master.x[(j,k)] for j in Parent.nRset)+
               sum(Parent.wZ[j]*Parent.Master.Y[(j,k)] for j in Parent.nZset))
    Parent.Master.DisjunctionBlock[k].BLOCK.cons.add(l_value - r_value >= 0)  #(82a) 
    
    for i in Parent.nLset:
        r_value = (Parent.s[i] - 
               sum(Parent.QR[(i,j)]*Parent.Master.xu[j] for j in Parent.mRset)-
               sum(Parent.QZ[(i,j)]*Parent.Master.yu[j] for j in Parent.mZset)-
               sum(Parent.PZ[(i,j)]*Parent.Master.Y[(j,k)] for j in Parent.nZset))
        l_value = sum(Parent.PR[(i,j)]*Parent.Master.x[(j,k)] for j in Parent.nRset)+ Parent.Master.wj[k]#(82b)
        
        Parent.Master.DisjunctionBlock[k].BLOCK.cons.add(r_value-l_value >= 0)
    
    #Dual Constraints
    
    for j in Parent.nRset:
        value= sum(Parent.Master.alpha[(i,k)]*Parent.PR[(i,j)] for i in Parent.nLset)
        value= value + sum(Parent.NR[(j,i)]*Parent.Master.beta[(i,k)] for i in Parent.zxset)  + Parent.Master.s[k,i]
        Parent.Master.DisjunctionBlock[k].BLOCK.cons.add(value >= Parent.wR[j])#(82c1)
    
    r_value=sum(Parent.Master.alpha[(j,k)] for j in Parent.nLset) - Parent.Master.tau[k]
    Parent.Master.DisjunctionBlock[k].BLOCK.cons.add(r_value == 0)
    
    r_value=sum(Parent.Master.beta[(j,k)]*Parent.Master.beta[(j,k)] for j in Parent.zxset)
    r_value=r_value+sum(Parent.Master.gamma[(j,k)]*Parent.Master.gamma[(j,k)] for j in Parent.zyset)
    Parent.Master.c_col.add(r_value <= Parent.Master.tau[k]*Parent.Master.tau[k])
    
    #Complementarity
    Parent.Master.DisjunctionBlock[k].BLOCK.comp1=ComplementarityList(rule=(complements(Parent.Master.alpha[(j,k)]>=0, 
                                  (Parent.s[j]-
                                   sum(Parent.QR[(j,i)]*Parent.Master.xu[i] for i in Parent.mRset)-
                                   sum(Parent.QZ[(j,i)]*Parent.Master.yu[i] for i in Parent.mZset)-
                                   sum(Parent.PR[(j,i)]*Parent.Master.x[(i,k)] for i in Parent.nRset)-
                                   sum(Parent.PZ[(j,i)]*Parent.Master.Y[(i,k)] for i in Parent.nZset) + Parent.Master.wj[k])>=0) for j in Parent.nLset)) #(82d)
    
    Parent.Master.DisjunctionBlock[k].BLOCK.comp2=ComplementarityList(rule=(complements(
            Parent.Master.x[(j,k)]>=0,
            Parent.s[(j,k)]>=0) for j in Parent.nRset)) #(82c2)
    
    #Conic Complementarity
    
    r_value=Parent.Master.tau[k]*Parent.Master.wj[k]-sum(Parent.Master.nxj[(j,k)]*Parent.Master.beta[(j,k)] for j in Parent.zxset) - sum(Parent.Master.nyj[(j,k)]*Parent.Master.gamma[(j,k)] for j in Parent.zyset)
    Parent.Master.c_col.add(r_value==0)
    for i in Parent.zxset:
        r_value=-Parent.Master.wj[k]*Parent.Master.beta[(i,k)]+Parent.Master.tau[k]*Parent.Master.nxj[(i,k)]
        Parent.Master.c_col.add(r_value==0)
    for i in Parent.zyset:
        r_value=-Parent.Master.wj[k]*Parent.Master.gamma[(i,k)]+Parent.Master.tau[k]*Parent.Master.nyj[(i,k)]
        Parent.Master.c_col.add(r_value==0)
        
    
    Parent.Master.DisjunctionBlock[k].c_disj=Disjunction(expr=[Parent.Master.DisjunctionBlock[k].LH, Parent.Master.DisjunctionBlock[k].BLOCK])
    #Parent.Master.pprint()
    
    return Parent

#Create parameters that will be updated on each iteration, initialized with coefficient vector of same size 
Parent.xu_star=Param(Parent.mRset,initialize=cR,default=0,mutable=True)
Parent.xl0_star=Param(Parent.nRset,initialize=wR,default=0,mutable=True)
Parent.yu_star=Param(Parent.mZset,initialize=cZ,default=0,mutable=True)
Parent.yl0_star=Param(Parent.nZset,initialize=wZ,default=0,mutable=True)
theta=0
Parent.theta=Param(initialize=theta,mutable=True)

Parent.xl_hat=Param(Parent.nRset,initialize=wR,within=NonNegativeReals,default=0,mutable=True)
Parent.yl_hat=Param(Parent.nZset,initialize=wZ,within=NonNegativeIntegers,default=0,mutable=True)

Parent.yl_star=Param(Parent.nZset, initialize=wZ,within=NonNegativeIntegers,default=0,mutable=True) 
Parent.xl_star=Param(Parent.nRset, initialize=wR,within=NonNegativeReals,default=0,mutable=True)
Parent.Theta_0=Param(initialize=0,mutable=True)
Parent.yl_arc=Param(Parent.nZset, initialize=wZ,within=NonNegativeIntegers,default=0,mutable=True)


''' Subproblem 1
    theta(xu*,yu*)=max wR*xl + wZ*yl
    s.t. PR*xl+PZ*yl <= s- QR*xu* - QZ*yu* (56)
'''
def sub1_obj(sub1):
    value=(sum(Parent.wR[j]*Parent.sub1.xl[j] for j in Parent.nRset)+
           sum(Parent.wZ[j]*Parent.sub1.yl[j] for j in Parent.nZset))
    return value

def sub1_c1(sub1,i):
    value=(sum(Parent.PR[(i,j)]*Parent.sub1.xl[j] for j in Parent.nRset)+
           sum(Parent.PZ[(i,j)]*Parent.sub1.yl[j] for j in Parent.nZset)+
           sum(Parent.QR[(i,j)]*Parent.xu_star[j] for j in Parent.mRset)+
           sum(Parent.QZ[(i,j)]*Parent.yu_star[j] for j in Parent.mZset))
    return Parent.zero <= Parent.s[i] -value -Parent.sub1.w

def sub1_c2(sub1,i):
    value=sum(Parent.NR[i,j]*Parent.sub1.xl[j] for j in Parent.nRset)
    return value==Parent.sub1.nx[i]

def sub1_c3(sub1,i):
    value=sum(Parent.NZ[i,j]*Parent.sub1.yl[j] for j in Parent.nZset)
    return value==Parent.sub1.ny[i]

def sub1_c4(sub1):
    value=sum(Parent.sub1.nx[j]*Parent.sub1.nx[j] for j in Parent.zxset )
    value=value+ sum(Parent.sub1.ny[j]*Parent.sub1.ny[j] for j in Parent.zyset)
    return value <= Parent.sub1.w*Parent.sub1.w


Parent.sub1=Block()   
Parent.sub1.xl=Var(Parent.nRset, within=NonNegativeReals)
Parent.sub1.yl=Var(Parent.nZset, within=NonNegativeIntegers)
Parent.sub1.w=Var(within=NonNegativeReals)
Parent.sub1.nx=Var(Parent.zxset,within=NonNegativeReals)
Parent.sub1.ny=Var(Parent.zyset, within=NonNegativeReals)
Parent.sub1.theta=Objective(rule=sub1_obj,sense=maximize)
Parent.sub1.c1=Constraint(Parent.nLset,rule=sub1_c1)
Parent.sub1.c2=Constraint(Parent.zxset, rule=sub1_c2)
Parent.sub1.c3=Constraint(Parent.zyset, rule=sub1_c3)
Parent.sub1.c4=Constraint(rule=sub1_c4)


''' Subproblem 2
    Theta_0(xu*,yu*)=min dR*xl + dZ*yl
    s.t. PR*xl+PZ*yl <= s- QR*xu* - QZ*yu* (56)
         BR*xl + BZ*yl <= r- AR*xu* - AZ*yu* (59)
         wR*xl + wZ*yl >= theta(xu*,yu*) (60)
'''
def sub2_obj(sub2):
    value=(sum(Parent.dR[i]*Parent.sub2.xl[i] for i in Parent.nRset)+
           sum(Parent.dZ[i]*Parent.sub2.yl[i] for i in Parent.nZset)) 
    return value

def sub2_c1(sub2,i):
    value=(sum(Parent.PR[(i,j)]*Parent.sub2.xl[j] for j in Parent.nRset)+
           sum(Parent.PZ[(i,j)]*Parent.sub2.yl[j] for j in Parent.nZset)+
           sum(Parent.QR[(i,j)]*Parent.xu_star[j] for j in Parent.mRset)+
           sum(Parent.QZ[(i,j)]*Parent.yu_star[j] for j in Parent.mZset))
    return Parent.zero <= Parent.s[i] -value

def sub2_c2(sub2,i):
    value=(sum(Parent.BR[(i,j)]*Parent.sub2.xl[j] for j in Parent.nRset)+
           sum(Parent.BZ[(i,j)]*Parent.sub2.yl[j] for j in Parent.nZset)+
           sum(Parent.AR[(i,j)]*Parent.xu_star[j] for j in Parent.mRset)+
           sum(Parent.AZ[(i,j)]*Parent.yu_star[j] for j in Parent.mZset))                                                                        
    return Parent.zero <= Parent.r[i] - value

def sub2_c3(sub2):
    value=(sum(Parent.wR[j]*Parent.sub2.xl[j] for j in Parent.nRset)+
           sum(Parent.wZ[j]*Parent.sub2.yl[j] for j in Parent.nZset))
    return value - Parent.theta >= Parent.zero


Parent.sub2=Block()
    
Parent.sub2.xl=Var(Parent.nRset, within=NonNegativeReals)
Parent.sub2.yl=Var(Parent.nZset, within=NonNegativeIntegers)
    
Parent.sub2.Theta_0=Objective(rule=sub2_obj,sense=minimize)

Parent.sub2.c1=Constraint(Parent.nLset,rule=sub2_c1) #(56)
Parent.sub2.c2=Constraint(Parent.mUset,rule=sub2_c2) #(59)
Parent.sub2.c3=Constraint(rule=sub2_c3)

def UBnew(Parent):
    UBnew=(sum(Parent.cR[j]*Parent.xu_star[j] for j in range(1,mR+1))+
        sum(Parent.cZ[j]*Parent.yu_star[j] for j in range(1,mZ+1))+
        Parent.Theta_0)
    return UBnew


#Iteration
while UB-LB > xi and k < maxit:
    #Step 1: Initialization (done)
    #Step 2: Solve the Master Problem 
    TransformationFactory('mpec.simple_disjunction').apply_to(Parent.Master)
    bigm.apply_to(Parent.Master) 
    #Parent.Master.Y.pprint()
    SCIP.solve(Parent.Master)
    if k > 1:
        print(f'sum of t ={sum(Parent.Master.t[(j,k)].value for j in Parent.nLset)}')
    for i in range(1,mR+1):
        Parent.xu_star[i]=Parent.Master.xu[i].value    
    for i in range(1,mZ+1):
        Parent.yu_star[i]=Parent.Master.yu[i].value    
    for i in range(1,nR+1):
        Parent.xl0_star[i]=Parent.Master.xl0[i].value
    for i in range(1,nZ+1):
        Parent.yl0_star[i]=Parent.Master.yl0[i].value

    LB=value(Parent.Master.Theta_star) 
    print(f'Iteration {k}: Master Obj={LB}')
    print(f'yu={Parent.Master.yu[1].value}')
    print(f'yl={Parent.Master.yl0[1].value}')
    #Step 3: Terminate?
    if UB-LB <= xi: #Output
        elapsed = time.time() - t
        flag=1 
        print(f'Optimal Solution Found in {k} iterations and {elapsed} seconds: Obj={UB}')
        break
    #Step 4: Solve first subproblem

    results1=GUROBI.solve(Parent.sub1) 
    
    if results1.solver.termination_condition !=TerminationCondition.optimal:
        raise RuntimeError("ERROR! ERROR! Subproblem 1: Could not find optimal solution")
    Parent.theta=value(Parent.sub1.theta)
    print(f'theta={value(Parent.theta)}')

    for i in range(1,nR+1):
        Parent.xl_hat[i]=Parent.sub1.xl[i].value 
    for i in range(1,nZ+1):
        Parent.yl_hat[i]=int(round(Parent.sub1.yl[i].value)) 

    
    #Step 5: Solve second subproblem
    results=GUROBI.solve(Parent.sub2)
    if results.solver.termination_condition==TerminationCondition.optimal: #If Optimal
        for i in range(1,nR+1):
            Parent.xl_star[i]=Parent.sub2.xl[i].value
        for i in range(1,nZ+1):
            Parent.yl_star[i]=int(round(Parent.sub2.yl[i].value))
            Parent.yl_arc[i]=int(round(Parent.sub2.yl[i].value))
        Parent.Theta_0=value(Parent.sub2.Theta_0)

 
        UB=min(UB,value(UBnew(Parent)))
        
    elif results.solver.termination_condition==TerminationCondition.infeasible or results.solver.termination_condition==TerminationCondition.infeasibleOrUnbounded: #If infeasible
        for i in range(1,nZ+1):
            Parent.yl_arc[i]=Parent.yl_hat[i]  
    else: 
         raise RuntimeError("ERROR! ERROR! Subproblem2 not infeasible or optimal solution not found: SOMETHING WENT VERY VERY WRONG") 
    
    #Step 6: Add new constraints
    k = k+1
    for i in range(1,nZ+1):
        Parent.Master.Y[(i,k)]=Parent.yl_arc[i] #Make sure yl_arc is int or else Master.Y rejects
    
    Master_add(Parent,k)
    #Step 7: Loop 


#Output Information regarding objective and time/iterations to convergence    
elapsed = time.time() - t
    
if k>= maxit:
    print('Maximum Iterations Reached')
elif k< maxit and flag !=1:
    print(f'Optimal Solution Found in {k-1} iterations and {elapsed} seconds: Obj={UB}')

