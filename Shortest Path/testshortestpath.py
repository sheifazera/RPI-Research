# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 11:53:50 2020

@author: sheifazera
"""
from pyomo.environ import *
import numpy as np
import random
from shortestpath import create_shortest_path, create_shortest_path_interdiction,create_shortest_path_robust_ellipsoid, return_path, return_path_robust
# %%
D={(0,1):22,(0,2):11,(1,3):69,(2,3):120,(1,0):22,(2,0):11,(3,1):69,(3,2):120}
N={0,1,2,3}
s=0
t=3

M=create_shortest_path(D,N,s,t)
opt=SolverFactory('gurobi')
results=opt.solve(M)
M.pprint()
path=return_path(M,D,N,s,t)
print(path)
#M.pprint()
# %% Interdiction

D={(0,1):(22,29),(0,2):(11,40),(1,3):(69,22),(2,3):(120,109)}
N={0,1,2,3}
s=0
t=3
B=1
M=create_shortest_path_interdiction(D, N, s, t, B)
opt=SolverFactory('gurobi')
results=opt.solve(M)
M.x.pprint()
print(f'd_t={M.d[t].value}')
path_interdicted=return_path_interdiction(M,D,N,s,t,opt)
print('Interdiction: optimal path')
print(path_interdicted)
#M.pprint()
# %%
#assume u is 3-dimensional so P is 4x3 dimensional also the edges go in the same order they do in the dictionary D...OR ELSE

D={(0,1):22,(0,2):11,(1,3):69,(2,3):120}
N={0,1,2,3}
s=0
t=3
udim=3
'''
#P_array=np.random.uniform(low=0, high=1, size=(4,udim))
#P={}
#j=0
#for (u,v) in D.keys():
#    for i in range(0,udim):
#        P[(u,v,i)]=P_array[j,i]
#    j=j+1
'''
P={(0,1,0):5,(0,1,1):40,(0,1,2):40,(0,2,0):15,(0,2,1):40,(0,2,2):5,(1,3,0):5,(1,3,1):0,(1,3,2):45,(2,3,0):45,(2,3,1):0,(2,3,2):5}
    
#P={(0,1,0):random.randrange(0,1),(0,1,1):random.randrange(0,1),(0,1,2):random.randrange(0,1),(0,1,3):random.randrange(0,1)}
# %%
M_ellipsoid=create_shortest_path_robust_ellipsoid(D, N, P, s, t, udim)
opt=SolverFactory('gurobi')
results_ellipsoid=opt.solve(M_ellipsoid)
path_E=return_path_robust(M_ellipsoid,D,N,P,s,t,udim)
print('Robust Shortest Path') 
print(path_E)
M_ellipsoid.pprint()
#M_ellipsoid.d.pprint()
#M_ellipsoid.u.pprint()

# %% Robustness test
udim=1
D={(0,1):22,(0,2):11,(1,3):69,(2,3):120}
s=0
t=3
N={0,1,2,3}
P={(0,1,0):(29),(0,2,0):(40),(1,3,0):(22),(2,3,0):(109)}

M_test=create_shortest_path_robust_ellipsoid(D,N,P,s,t,udim)
test=opt.solve(M_test)
print('Test: Robust Shortest Path') 
M_test.pprint()
#M_test.d.pprint()



# %% Robust Interdiction
udim=3
s=0
t=3
P={(0,1,0):5,(0,1,1):40,(0,1,2):40,(0,2,0):15,(0,2,1):40,(0,2,2):5,(1,3,0):5,(1,3,1):0,(1,3,2):45,(2,3,0):45,(2,3,1):0,(2,3,2):5}
D={(0,1):(22,29),(0,2):(11,40),(1,3):(69,22),(2,3):(120,109)}
N={0,1,2,3}
B=1

M_RI=create_shortest_path_robust_ellipsoid_certain_interdiction(D,N,P,s,t,udim,B)
opt=SolverFactory('gurobi')
RI=opt.solve(M_RI,tee=True)
path_RI=return_path_interdiction_robust(M_RI,D,N,P,s,t,udim,opt)
print(f'd_t={M_RI.d[t].value}')
print('Robust Shortest Path with Interdiction')
print(path_RI)
M_RI.x.pprint()
M_RI.d.pprint()
M_RI.u.pprint()

#print('Robust SP: Certain Interdiction Lengths')
#M_RI.pprint()


# %% Robust Interdiction Uncertain Additions
udim=3
s=0
t=3
D={(0,1):(22,29),(0,2):(11,40),(1,3):(69,22),(2,3):(120,109)}
P={(0,1,0):5,(0,1,1):40,(0,1,2):40,(0,2,0):15,(0,2,1):40,(0,2,2):5,(1,3,0):5,(1,3,1):0,(1,3,2):45,(2,3,0):45,(2,3,1):0,(2,3,2):5}
N={0,1,2,3}
B=1
vdim=2
R={(0,1,0):10,(0,1,1):1,(0,2,0):10,(0,2,1):1,(1,3,0):1,(1,3,1):15,(2,3,0):1,(2,3,1):1} #norm(R)=15.9146
M_RIU=create_shortest_path_robust_ellipsoid_uncertain_interdiction(D,N,P,R,s,t,udim,vdim,B)
opt=SolverFactory('gurobi')
RIU=opt.solve(M_RIU)
path_RIU=return_path_interdiction_robust_uncertain(M_RIU,D,N,P,R,s,t,udim,vdim,opt)
print(f'd_t={M_RIU.d[t].value}')
print('Robust Shortest Path with Interdiction and Uncertain Interdiction Additions')
print(path_RIU)
M_RIU.x.pprint()
M_RIU.d.pprint()
M_RIU.u.pprint()


# %% 
s=0
t=3
opt=SolverFactory('gurobi')
N={0,1,2,3}
B=2
D={(0,1):(2,3) , (0,2): (3,1), (1,3): (3,4), (2,3): (3,2)}

vdim=2
R={(0,1,0): 1 , (0,1,1): 0 , (0,2,0): 0, (0,2,1): 2, (1,3,0): 3, (1,3 ,1): 0, (2,3,0): 0, (2,3 ,1): 3}

M=create_asymmetric_uncertainty_shortest_path_interdiction(D,N,R,s,t,vdim,B)
results=opt.solve(M,tee=True)
print('Asymmetric Uncertainty Interdiction Results')
print('Interdiction Variable')
M.x.pprint()
print('Evader Path Variable')
M.w.pprint()
M.z.pprint()
print('Uncertainty Variable')
M.t.pprint()
path=return_path_asymmetric(M,D,N,s,t)
print(path)

# %%
s=0
t=3
N={0,1,2,3}
B=4
Prob={(0,1):(0.8,0.6),(0,2):(0.9,0.6), (1,3):(0.5,0.4), (2,3): (0.3,0.1)}
D={}
for (i,j) in Prob.keys():
    (p,q)=Prob[(i,j)]
    D[(i,j)]=(-np.log(p),-np.log(q)+np.log(p))

vdim=2
R={(0,1,0):0.29, (0,1,1): 0, (0,2,0):0, (0,2,1):.4, (1,3,0):0.22, (1,3,1):0, (2,3,0):0, (2,3,1):1.09}
M=create_asymmetric_uncertainty_shortest_path_interdiction(D,N,R,s,t,vdim,B)
opt=SolverFactory('gurobi')
results=opt.solve(M)
print('Asymmetric Uncertainty Interdiction Results')
print('Interdiction Variable')
M.x.pprint()
print('Evader Path Variable')
M.w.pprint()
M.z.pprint()
print('Uncertainty Variable')
M.t.pprint()
path=return_path_asymmetric(M,D,N,s,t)
print(path)
print(f'Objective={value(M.Obj)}')
print(f'd_t={M.d[t].value}')


