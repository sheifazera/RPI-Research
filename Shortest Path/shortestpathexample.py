# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 10:25:33 2020

@author: spunlag
"""
from pyomo.environ import *


D={(0,1):22,(0,2):11,(1,3):69,(2,3):120}
s=0
t=3

M=create_shortest_path(D,s,t)