#Sioux Falls Network
import numpy as np
from pyomo.environ import *
from shortestpath import create_shortest_path, create_shortest_path_robust_ellipsoid
matrix = np.loadtxt('siouxfalls.txt', usecols=(0,1,3))
#print(matrix) 

D = {(int(a),int(b)) : c for a,b,c in matrix}
N=range(1,25)
s=11
t=20
#print(D)

M=create_shortest_path(D,N,s,t)
opt=SolverFactory('gurobi')
results=opt.solve(M)
M.pprint()



ME=create_shortest_path_robust_ellipsoid(D,N,P,s,t)
results_robust_ellipsoid=opt.solve(ME)
ME.pprint()
