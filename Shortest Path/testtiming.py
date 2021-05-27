# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 09:39:05 2021

@author: sheifazera
"""
from pyomo.environ import *
M=ConcreteModel()
M.x=Var()
M.c1=Constraint(expr=M.x+1 <= 4)
M.obj=Objective(expr=M.x, sense=maximize)
opt=SolverFactory('gurobi')
results=opt.solve(M,tee=True)

