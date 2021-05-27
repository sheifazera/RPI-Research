

from pao.mpr import *
 

 

def create():
    M = LinearMultilevelProblem()
 

    U = M.add_upper(nxZ=1)
    L = U.add_lower(nxZ=1)
 

    # Variables
    U.x.lower_bounds = [0.0]
    U.x.upper_bounds = [10.0]
    L.x.lower_bounds = [0.0]
    L.x.upper_bounds = [5.0]
 

    # Objectives
    U.c[U] = [-1]
    U.c[L] = [-7]
 

    L.c[L] = [1]
 

    # Constraints
    U.A[U] = [[-3], [1]]
    U.A[L] = [[2], [2]]
    U.b = [12, 20]
 

    L.A[U] = [[2], [-2]]
    L.A[L] = [[-1], [4]]
    L.b = [7, 16]
 

    return M
 

if __name__ == "__main__":          #pragma: no cover
    M = create()
    M.check()
    #M.print()
    opt = Solver('pao.mpr.PCCG')
    opt.solve(M, quiet=False, mip_solver='gurobi')
    print(M.U.x.values)
    print(M.U.LL.x.values)
    
 # She'ifa's Check   
    
from pyomo.environ import * 

M_check=ConcreteModel()
M.yU=Var(within=Integers)
M.Obj=Objective(expr=-M.yU - 35, sense=minimize )
M.c1=

