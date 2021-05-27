# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 10:26:08 2021

@author: sheifazera
"""
import itertools
import random


N=list(range(1,1501))
all_arcs=list(itertools.permutations(N,2))

seeds={631996,5111994,7161951,7171956,12152015,2122021,2142021,2162021,2182021,2282021}
for A in seeds:
    random.seed(a=A)
    arcs_1500=random.choices(all_arcs,k=7500)
    Prob={}
    for (i,j) in arcs_1500:
        Prob[(i,j)]=random.random()
        
        
        




