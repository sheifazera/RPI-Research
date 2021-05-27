# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 10:28:36 2021

@author: sheifazera
"""

import networkx as nx
from pyomo.environ import *
import numpy as np
from shortestpath_networkx import create_shortest_path_nx, create_shortest_path_interdiction_nx, create_asymmetric_uncertainty_shortest_path_interdiction_nx


G=nx.DiGraph()
G.add_edge(1, 2, length=0.8 )
G.add_edge(1, 3, length=0.6)
G.add_edge(2, 4, length=0.1)
G.add_edge(3, 4, length=0.8)


M=create_shortest_path_nx(G,1,4)
opt=SolverFactory('gurobi')
opt.solve(M)

#%%
G=nx.DiGraph()
G.add_edge(1,2, length=0.8, interdicted_length=0.4)
G.add_edge(1,3, length=0.6, interdicted_length=0.2)
G.add_edge(2,4, length=0.1, interdicted_length=0.3)
G.add_edge(3,4, length=0.8, interdicted_length=0.2)

M=create_shortest_path_interdiction_nx(G,1,4,1)
opt.solve(M)


#%% 
G=nx.DiGraph()
Prob={(0,1):(0.8,0.6),(0,2):(0.9,0.6), (1,3):(0.5,0.4), (2,3): (0.3,0.1)}
D={}
for (i,j) in Prob.keys():
    (p,q)=Prob[(i,j)]
    G.add_edge(i,j,length=-np.log(p), interdicted_length=-np.log(q)+np.log(p))
R={(0,1,0):0.29, (0,1,1): 0, (0,2,0):0, (0,2,1):.4, (1,3,0):0.22, (1,3,1):0, (2,3,0):0, (2,3,1):1.09}
vdim=2
B=1
M=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R,0,3,vdim,B)
opt.solve(M)

#%%

