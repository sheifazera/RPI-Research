# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 09:35:59 2021

@author: sheifazera
"""


import numpy as np
import itertools as it
import matplotlib.pyplot as plt
import random
import time
from pyomo.environ import *
from prettytable import PrettyTable, ALL
import networkx as nx
import pickle

#%%

Nodes=[10,100,500,1000]
edges=np.empty([20,4])

for loop in range(1,21):
    l=0
    for N in Nodes:
        connected=False
        while connected is False:
            G=nx.random_geometric_graph(N,1.5/sqrt(N), dim=2, p=2)
            connected=nx.is_connected(G)
        '''
        plt.subplot(121)
        nx.draw(G, with_labels=False, font_weight='bold', node_size=25)
        '''
        edges[loop-1][l]=len(list(G.edges))
        l=l+1
    print(f'Finished loop {loop}') 
#%%
set_label=np.array([10,100,500,1000])
plt.figure()
plt.title('Number of Edges')
plt.xlabel('Nodes')
plt.ylabel('Edges')
plt.boxplot(edges,labels=set_label)


#%%
