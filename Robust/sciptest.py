# -*- coding: utf-8 -*-
"""
Created on Tue May  5 12:36:24 2020

@author: sheifazera
"""

import time
import sys 

from pyomo.environ import *
from pyomo.gdp import *
from pyomo.mpec import *


solvername='SCIPAMPL'

solverpath_folder='C:\\Users\sheifazera\Documents\RPI\Research\SCIP' 

solverpath_exe='C:\\Users\sheifazera\Documents\RPI\Research\SCIP\scipampl'

sys.path.append(solverpath_folder)

SCIP=SolverFactory(solvername,executable=solverpath_exe)

m=ConcreteModel()

m.x=Var(within=NonNegativeReals)
m.y=Var(within=NonNegativeReals)
m.z=Var(within=NonNegativeReals)

def Obj(m):
    value = m.x + m.y + m.z
    return value
m.Obj=Objective(rule=Obj,sense=minimize)

def NonlinearConstraint(m):
    value=m.x*m.y*m.y
    return value >= 1

m.NonlinearConstraint=Constraint(rule=NonlinearConstraint)
m.pprint()
results=SCIP.solve(m,tee=True)

print(f'x={m.x.value},y={m.y.value},z={m.z.value}')