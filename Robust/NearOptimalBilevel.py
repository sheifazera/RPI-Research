# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 10:08:12 2020

@author: She'ifa
"""

import time

from pyomo.environ import *
from pyomo.gdp import *
from pyomo.mpec import *

infinity = float('inf')
opt=SolverFactory("gurobi")
bigm = TransformationFactory('gdp.bigm')


from Parameters import 

#Solve High Point Relaxation

#Solve Optimistic Bilevel

#Dual Feasibility (Problem 5.3)

#Assemble vertices

#Near-optimal Bilevel (Problem 5.6 and 5.7)