import numpy as np
from pyomo.environ import *
import networkx as nx
#from shortestpath import create_shortest_path, create_shortest_path_robust_ellipsoid, return_path, return_path_robust #, return_path_robust
from shortestpath_networkx import return_paths, create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx, return_interdicted_arcs, return_paths_multiple_ST, return_interdicted_arcs, create_asymmetric_uncertainty_shortest_path_interdiction_nx
#%%
matrix = np.loadtxt('siouxfalls.txt', usecols=(0,1,3))

D = {(int(a),int(b)) : c for a,b,c in matrix}
N=set(range(1,25))
s=3
t=20

Prob = {(int(a),int(b)) : c for a,b,c in matrix}
for (i,j) in Prob.keys():
    if Prob[(i,j)] ==10:
        Prob[(i,j)] = (0.1,0.05) #p_ij, q_ij 
    elif Prob[(i,j)] ==9 :
        Prob[(i,j)] = (0.2, 0.1)
    elif Prob[(i,j)] ==8 :  
        Prob[(i,j)] = (0.3,0.15)
    elif Prob[(i,j)] ==7 :  
        Prob[(i,j)] = (0.4,0.2)
    elif Prob[(i,j)] ==6 :  
        Prob[(i,j)] = (0.5, 0.25)
    elif Prob[(i,j)] ==5 :  
        Prob[(i,j)] = (0.6, 0.3)
    elif Prob[(i,j)] ==4 :  
        Prob[(i,j)] = (0.7, 0.35)
    elif Prob[(i,j)] ==3 :  
        Prob[(i,j)] = (0.8, 0.4)
    elif Prob[(i,j)] ==2 :  
        Prob[(i,j)] = (0.9, 0.45)    
    else:
        raise Exception('Original Sioux Falls value out of range 2-10')
G=nx.DiGraph()
for (i,j) in Prob.keys():
    (p,q)=Prob[(i,j)]
    G.add_edge(i,j,length=-np.log(p),interdicted_length=-np.log(q)+np.log(p))
    #r=0.25*q
    #G.add_edge(i,j,length=-np.log(p),interdicted_length=-np.log(q)+np.log(p), robust=-np.log(r)+np.log(q))
vdim=len(set(G.edges()))
R={}

for (i,j) in set(G.edges):
    r=G[i][j]['interdicted_length']
    k=0
    for (u,v) in set(G.edges):
        if i==u and j==v:
            R[(i,j,k)]=r
        else:
            R[(i,j,k)]=0
        k=k+1    

N_external={1,2,3,7,12,13,18,20,21,24}
I=set()
for s in N_external:
    I.add((s,10))

B=5
opt=SolverFactory('gurobi_direct')


def get_IR(M,G,R,s,t,vdim,B):
    #Solve the nominal model
    R0={}
    vdim0=0
    M_nom=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R0, s,t,vdim0,B)
    opt.solve(M_nom)
    #(paths,length)=return_paths(M_nom,G,R0,s,t)
    #print('Nominal SPI: Path')
    #print(paths)
    #print(length)
    #Solve a model that has the fixed interdictions from the nominal model but adds in robustness
    M_nom_rob=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R, s,t,vdim,B)
    M_nom_rob.mods=ConstraintList()
    #print('Nominal SPI: Interdicted Arcs')
    for (i,j) in G.edges:
        if value(M_nom.x[(i,j)]) >= 0.9:
            M_nom_rob.mods.add(M_nom_rob.x[(i,j)]==1)
            #print(f'({i},{j})')
        else:
            M_nom_rob.mods.add(M_nom_rob.x[(i,j)]==0)
    opt.solve(M_nom_rob)
    #M_nom_rob.x.pprint()
    #(paths,length)=return_paths(M_nom_rob,G,R,s,t)
    #print('Robust Path Given Nominal Interdictions')
    #print(paths)
    #print(length)
    #Make comparisons
    #print(value(M.Obj))
    #print(value(M_nom.Obj))
    #print(value(M_nom_rob.Obj))
    Regret=100*(exp(-value(M.Obj))-exp(-value(M_nom.Obj)))/exp(-value(M.Obj))
    Improvement=100*(exp(-value(M_nom_rob.Obj))-exp(-value(M.Obj)))/exp(-value(M_nom_rob.Obj))

    return (Improvement,Regret)




def get_IR_multiple_ST(M,G,R,I,S,T,vdim,B):
    R0={}
    vdim0=0
    M_nom=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G,R0, I ,vdim0,B)
    opt.solve(M_nom)
    M_nom_rob=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G, R, I, vdim, B)
    M_nom_rob.mods=ConstraintList()
    for (i,j) in G.edges:
        if value(M_nom.x[(i,j)]) >= 0.9:
            M_nom_rob.mods.add(M_nom_rob.x[(i,j)]==1)
            #print(f'({i},{j})')
        else:
            M_nom_rob.mods.add(M_nom_rob.x[(i,j)]==0)
    opt.solve(M_nom_rob)
    (_,length_nom_rob)=return_paths_multiple_ST(M_nom_rob, G, R, S, T)
    (_,length)=return_paths_multiple_ST(M, G, R, S, T)
    (_,length_nom)=return_paths_multiple_ST(M_nom, G,R,S,T)
    Regret_Avoided=100*(exp(-length_nom_rob)-exp(-length))/exp(-length_nom_rob)
    Objective_Paid=100*(exp(-length)-exp(-length_nom))/exp(-length)
    
    return(Regret_Avoided,Objective_Paid)
    

#%% Super Source

G_SS=G
for i in N_external:
    #G_SS.add_edge(i,100, length=-np.log(1),interdicted_length=-np.log(1)+np.log(1))
    G_SS.add_edge(100,i, length=-np.log(1),interdicted_length=-np.log(1)+np.log(1))
 
vdim_SS=len(set(G_SS.edges()))
R_SS={}

for (i,j) in set(G_SS.edges):
    r=G[i][j]['interdicted_length']
    k=0
    for (u,v) in set(G_SS.edges):
        if i==u and j==v:
            R_SS[(i,j,k)]=r
        else:
            R_SS[(i,j,k)]=0
        k=k+1  
S=100     
T=10        
M_SSSPIAU=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G_SS, R_SS,S,T, vdim_SS, B)
opt.solve(M_SSSPIAU)
interdicted=return_interdicted_arcs(M_SSSPIAU,G_SS)
(paths,lengths)=return_paths(M_SSSPIAU,G_SS,R_SS,S,T)
(Regret_Avoided,Objective_Paid)=get_IR(M_SSSPIAU,G_SS,R_SS,S,T,vdim_SS,B)
print('SS-SPIAU')
print('Interdict')
print(interdicted)
print('Path')
print(paths)
print('Adjusted Length')
print(np.exp(-lengths))
print('Evader Length')
print(np.exp(-M_SSSPIAU.d[T].value))
print('Regret Avoided')
print(Regret_Avoided)

R_SS={}
vdim_SS=0
M_SSSPI=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G_SS, R_SS, S,T, vdim_SS, B)
opt.solve(M_SSSPI)
interdicted=return_interdicted_arcs(M_SSSPI,G_SS)
(paths,lengths)=return_paths(M_SSSPI,G_SS,R_SS,S,T)
(Regret_Avoided,Objective_Paid)=get_IR(M_SSSPI,G_SS,R_SS,S,T,vdim_SS,B)
print('SS-SPI')
print('Interdict')
print(interdicted)
print('Path')
print(paths)
print('Adjusted Length')
print(np.exp(-lengths))
print('Evader Length')
print(np.exp(-M_SSSPI.d[T].value))
print('Regret Avoided')
print(Regret_Avoided)
#%% Worst-case (S,T)

#random.seed(a=631996)

R_WC={}
vdim_WC=len(set(G.edges))

for (i,j) in set(G.edges):
    r=G[i][j]['interdicted_length']
    k=0
    for (u,v) in set(G.edges):
        if i==u and j==v:
            R_WC[(i,j,k)]=r
        else:
            R_WC[(i,j,k)]=0
        k=k+1 
        
M_WCSPIAU=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G, R_WC, I, vdim_WC, B)
opt.solve(M_WCSPIAU)
interdicted=return_interdicted_arcs(M_WCSPIAU,G)
print('WC-SPIAU')
print('Interdict')
print(interdicted)
for (S,T) in I:
    (paths,lengths)=return_paths_multiple_ST(M_WCSPIAU,G,R_WC,S,T)
    (Regret_Avoided,Objective_Paid)=get_IR_multiple_ST(M_WCSPIAU,G,R_WC,I, S,T,vdim_WC,B)
    print(f'({S},{T})')
    print('Path')
    print(paths)
    print('Adjusted Length')
    print(np.exp(-lengths))
    print('Evader Length')
    print(np.exp(-M_WCSPIAU.STBlocks[(S,T)].d[T].value))
    print('Regret Avoided')
    print(Regret_Avoided)

#%%
R_WC={}
vdim_WC=0
M_WCSPI=create_asymmetric_uncertainty_shortest_path_interdiction_multiple_ST_nx(G, R_WC, I, vdim_WC, B)
opt.solve(M_WCSPI)
interdicted=return_interdicted_arcs(M_WCSPI,G)
print('WC-SPI')
print('Interdict')
print(interdicted)
for (S,T) in I:
    (paths,lengths)=return_paths_multiple_ST(M_WCSPI,G,R_WC,S,T)
    (Regret_Avoided,Objective_Paid)=get_IR_multiple_ST(M_WCSPI,G,R_WC,I, S,T,vdim_WC,B)
    print(f'({S},{T})')
    print('Path')
    print(paths)
    print('Adjusted Length')
    print(np.exp(-lengths))
    print('Evader Length')
    print(np.exp(-M_WCSPI.STBlocks[(S,T)].d[T].value))
    print('Regret Avoided')
    print(Regret_Avoided)


