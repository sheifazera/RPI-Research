# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 10:29:21 2021

@author: sheifazera
"""

import numpy as np
import itertools as it
import matplotlib.pyplot as plt
import random
import time
from pyomo.environ import *
from prettytable import PrettyTable, ALL
from shortestpath_networkx import create_shortest_path_nx, create_shortest_path_interdiction_nx, create_asymmetric_uncertainty_shortest_path_interdiction_nx, return_paths, similar_length, create_shortest_path_interdiction_nx
import networkx as nx
import pickle
#%%
matrix = np.loadtxt('eastern_mass.txt', usecols=(0,1,3))
Prob = {(int(a),int(b)) : c for a,b,c in matrix}

A=max(Prob.values())
B=min(Prob.values())
a=1
b=0.1

for (i,j) in Prob.keys():
    Prob[(i,j)]=(a-b)*(Prob[(i,j)]-B)/(A-B) + b
    
G=nx.DiGraph()
for (i,j) in Prob.keys():
    p=Prob[(i,j)]
    q=0.5*p
    G.add_edge(i,j,length=-np.log(p),interdicted_length=-np.log(q)+np.log(p))

vdim=len(set(G.edges()))
R={}

for (i,j) in set(G.edges):
    l=G[i][j]['length']
    r=G[i][j]['interdicted_length']
    k=0
    for (u,v) in set(G.edges):
        if i==u and j==v:
            R[(i,j,k)]=r
        else:
            R[(i,j,k)]=0
        k=k+1    
matrix=None
Prob=None  

B=25 #approximately 10 percent of the arcs can be interdicted

I=set(it.combinations(list(G.nodes),2))
data={}
cutoff=1.1
#data (S,T): {Djikstra solve time, SP solve time, SPI solve time, Robust solve time, path shortest path length, number of nodes in shortest path, number of similar paths}     
#%% 
table=PrettyTable()
table.hrules=ALL
table.field_names=["(S,T)", "Solve Time (s)", "Path", "Shortest Path Length", "Number of Nodes in Shortest Path"]
loop=0
   
opt=SolverFactory('gurobi')
opt.options['TimeLimit']=600
for (S,T) in I:
    print(f'Percentage Complete: {100*loop/2700}')
    start=time.time()
    paths=nx.shortest_path(G, S, T, weight='length')
    end=time.time()
    djikstras_time=end-start
    M_SP=create_shortest_path_nx(G,S,T)
    start=time.time()
    opt.solve(M_SP)
    end=time.time()
    SP_time=end-start
    M_interdiction=create_shortest_path_interdiction_nx(G, S,T, B)
    start=time.time()
    opt.solve(M_interdiction)
    end=time.time()
    interdiction_time=end-start
    M=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R,S,T,vdim,B)   
    start=time.time()
    opt.solve(M)
    end=time.time()
    (paths, lengths)=return_paths(M,G,R,S,T) #generator
    #To do once we determine how shortest path returns paths
    #if len(paths) > 1:
    '''
    path_cell=[]
    length_cell=[]
    number_cell=[]
    for path in paths:
        path_str=[]
        for i in path:
            path_str.append(str(i))
        path_str="-".join(path_str)
        path_cell.append(path_str)
        length_cell.append(str(lengths))
        number_cell.append(str(len(path)))
    path_cell="\n".join(path_cell)    
    length_cell="\n".join(length_cell)    
    number_cell="\n".join(number_cell)
    
    
    table.add_row([(S,T), djikstras_time, SP_time, interdiction_time, end-start, path_cell,length_cell, number_cell])
    
    #else:
        #len
        #table.add_row([(S,T), end-start, paths[0], lengths, len(paths[0])])
        
        '''
    (similar_paths, no_of_similar_paths)=similar_length(M,G,R,S,T, cutoff)
    data[(S,T)]=(djikstras_time, SP_time, interdiction_time, end-start,  paths[0], lengths, len(paths[0]), no_of_similar_paths)
    loop=loop+1     

filename='ST_data_with_similar'
outfile = open(filename,'wb')
pickle.dump(data,outfile)
outfile.close()
'''
table_txt=table.get_string()
with open('ST_data_table.txt','a') as file:
    file.write(table_txt)
'''
#%% Figures (Matplotlib Scatter plots?)
pickle_off = open("ST_data_with_similar","rb")
data = pickle.load(pickle_off)


k=0
time_array=np.empty([2701,1])
length_array=np.empty([2701,1])
number_array=np.empty([2701,1])

time_0_1=[]
time_1_2=[]
time_2_3=[]
time_3_4=[]
time_4_5=[]
time_5_6=[]
time_6_7=[]
time_7_8=[]
time_8_9=[]
time_9_10=[]
time_10_11=[]
time_11_12=[]
time_12_13=[]
time_13_14=[]
time_14_15=[]
time_15_16=[]
time_16_17=[]

djikstra_array=np.empty([2701,1])
SP_array=np.empty([2701,1])
SPI_array=np.empty([2701,1])
similar_array=np.empty([2701,1])

for (S,T) in data.keys():
    (djikstra, SP, SPI, t , path , l , n, no_of_similar)=data[(S,T)]
    
    djikstra_array[k,0]=djikstra
    SP_array[k,0]=SP
    SPI_array[k,0]= SPI
    time_array[k,0]=t
    length_array[k,0] =l
    number_array[k,0] =n
    similar_array[k,0] = no_of_similar
    
    if l <1 and l >=0:
        time_0_1.append(t)
    elif l<2 and l>=1:
        time_1_2.append(t)
    elif l<3 and l>=2:
        time_2_3.append(t)
    elif l<4 and l>=3:
        time_3_4.append(t)
    elif l<5 and l>=4:
        time_4_5.append(t)
    elif l<6 and l>=5:
        time_5_6.append(t)
    elif l<7 and l>=6:
        time_6_7.append(t)
    elif l<8 and l>=7:
        time_7_8.append(t)
    elif l<9 and l>=8:
        time_8_9.append(t)
    elif l<10 and l>=9:
        time_9_10.append(t)
    elif l<11 and l>=10:
        time_10_11.append(t)
    elif l<12 and l>=11:
        time_11_12.append(t)
    elif l<13 and l>=12:
        time_12_13.append(t)
    elif l<14 and l>=13:
        time_13_14.append(t)
    elif l<15 and l>=14:
        time_14_15.append(t)
    elif l<16 and l>=15:
        time_15_16.append(t)
    elif l<17 and l>=16:
        time_16_17.append(t)
    
    k=k+1
    

fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(length_array, time_array)
ax.set_xlabel('Shortest Path Length')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Length vs Time')

fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(number_array, time_array)
ax.set_xlabel('Number of Nodes in the Shortest Path')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Nodes vs Time')


fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(similar_array, time_array)
ax.set_xlabel('Number Similar Length Paths')
ax.set_ylabel('Solve Time (s)')
ax.set_title('Number of Similar Length Paths (10%)')

#%%
fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.hist(length_array, bins=19, range=(0,18))
ax.set_xlabel('Shortest Path Length')
ax.set_ylabel('Occurences')

fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.hist(length_array, bins=17, range=(0,16))
ax.set_xlabel('Number of Nodes in Shortest Path')
ax.set_ylabel('Occurences')

fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.hist(length_array, bins=31, range=(0,30))
ax.set_xlabel('Number of Similar Length Paths')
ax.set_ylabel('Occurences')



#%%
binned_data=[np.asarray(time_0_1), np.asarray(time_1_2),np.asarray(time_2_3),np.asarray(time_3_4),np.asarray(time_4_5),np.asarray(time_5_6), np.asarray(time_6_7), np.asarray(time_7_8),np.asarray(time_8_9), np.asarray(time_9_10), np.asarray(time_10_11), np.asarray(time_11_12), np.asarray(time_12_13), np.asarray(time_13_14), np.asarray(time_14_15), np.asarray(time_15_16), np.asarray(time_16_17)]


fig1, ax1 = plt.subplots()
ax1.set_title('Solve time by length of shortest path')
ax1.boxplot(binned_data, labels=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])
ax1.set_xlabel('Shortest Path Length')
ax1.set_ylabel('Solve Time (s)')

fig1, ax1 = plt.subplots()
ax1.set_title('Solve time by length of shortest path')
ax1.boxplot(binned_data, labels=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17], showfliers=False)
ax1.set_xlabel('Shortest Path Length')
ax1.set_ylabel('Solve Time (s)')


#%% 

fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(length_array, similar_array)
ax.set_ylabel('Number Similar Length Paths (10%)')
ax.set_xlabel('Shortest Path Length')
ax.set_title('Number of Similar Length Paths and Shortest Path Length')


#%%
fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(number_array, similar_array)
ax.set_ylabel('Number Similar Length Paths (10%)')
ax.set_xlabel('Nodes in Shortest Path')
ax.set_title('Number of Similar Length Paths and Nodes in Shortest Path')

#%%
Budgets={5,7,10,12,15,17,20}
Budgets={25}
opt=SolverFactory('gurobi')
opt.options['TimeLimit']=600

S=4
T=36
'''
paths=nx.shortest_path(G, S, T, weight='length')

M_interdiction=create_shortest_path_interdiction_nx(G, S,T, B)
opt.solve(M_interdiction, tee=True)
M_interdiction.x.pprint()
'''

for B in Budgets:
    M=create_asymmetric_uncertainty_shortest_path_interdiction_nx(G,R,S,T,vdim,B)
    start=time.time()
    opt.solve(M)
    end=time.time()
    print(f'Time to solve at budget {B} is {end-start}')
    (paths, length) =return_paths(M,G,R,S,T)
    for path in paths:
        print(path)
        print(length)
    for (i,j) in G.edges:
        if M.x[(i,j)].value >= 0.9:
            print(f'({i},{j})')



#%% Heat Map for SP Length

from matplotlib.colors import LogNorm
x_min = np.min(length_array[:,0]) 
x_max = np.max(length_array[:,0]) 
  
y_min = np.min(similar_array[:,0]) 
y_max = np.max(similar_array[:,0]) 
  
x_bins = np.linspace(x_min, x_max, 19) 
y_bins = np.linspace(y_min, y_max, 31)

fig, ax = plt.subplots(figsize =(10, 7)) 
plt.hist2d(length_array[:,0], similar_array[:,0], bins =[x_bins, y_bins],
                norm=LogNorm())
ax.set_ylabel('Number Similar Length Paths (10%)')
ax.set_xlabel('Shortest Path Length')
ax.set_title('Number of Similar Length Paths and Shortest Path Length')
plt.colorbar() 
#%% Heat Map for Nodes in SP


x_min = np.min(number_array[:,0]) 
x_max = np.max(number_array[:,0]) 
  
y_min = np.min(similar_array[:,0]) 
y_max = np.max(similar_array[:,0]) 
  
x_bins = np.linspace(x_min, x_max, 15) 
y_bins = np.linspace(y_min, y_max, 31)

fig, ax = plt.subplots(figsize =(10, 7)) 
plt.hist2d(number_array[:,0], similar_array[:,0], bins =[x_bins, y_bins], norm=LogNorm())
                
ax.set_ylabel('Number Similar Length Paths (10%)')
ax.set_xlabel('Nodes in Shortest Path')
ax.set_title('Number of Similar Length Paths and Nodes in Shortest Path')
plt.colorbar()




#%% Multivariate Regression
import pandas as pd

from sklearn import linear_model

from scipy.stats import pearsonr



pickle_off = open("ST_data_with_similar","rb")
data = pickle.load(pickle_off)
k=0
data_array=np.zeros(shape=(2701,4))
for (S,T) in data.keys():
    (_,_,_,t,_,l,n,s)=data[(S,T)]
    data_array[k]=np.array([t,l,n,s])
    k=k+1
    
#data_array = np.array([list(v) for v in data.values()])
X=pd.DataFrame(data_array[:,1:4], columns=['Path_Length','Nodes_in_Path', 'Similar_Paths'])
y=pd.DataFrame(data_array[:,0], columns=['Solve_Time'])

#%%                                           
regr = linear_model.LinearRegression()
model=regr.fit(X, y)
response = model.predict(X)
r2=model.score(X,y)


print(regr.coef_)  
print(r2) 

X.corr


# Path_Length
X_path=pd.DataFrame(data_array[:, 1], columns=['Path_Length'])

regr_path=linear_model.LinearRegression()
model=regr_path.fit(X_path,y)
response= model.predict(X_path)
r2=model.score(X_path,y)   
print('Path Length vs Solve Time')
print(regr_path.coef_)
print(r2) 

plt.style.use('default')
plt.style.use('ggplot')

fig, ax = plt.subplots(figsize=(8, 4))

ax.plot(X_path, response, color='k', label='Regression model')
ax.scatter(X_path, y, edgecolor='k', facecolor='grey', alpha=0.7, label='Sample data')
ax.set_ylabel('Solve Time (s)', fontsize=14)
ax.set_xlabel('Path Length', fontsize=14)
ax.legend(facecolor='white', fontsize=11)
ax.set_title('$R^2= %.2f$' % r2, fontsize=18)

fig.tight_layout()

#Similar Paths
X_sim=pd.DataFrame(data_array[:,3], columns=['Similar_Paths'])

regr_sim=linear_model.LinearRegression()
model=regr_path.fit(X_sim,y)
response= model.predict(X_sim)
r2=model.score(X_sim,y)   
print('Number of Similar Paths vs Solve Time')
print(regr_path.coef_)
print(r2)

fig, ax = plt.subplots(figsize=(8, 4))

ax.plot(X_sim, response, color='k', label='Regression model')
ax.scatter(X_sim, y, edgecolor='k', facecolor='grey', alpha=0.7, label='Sample data')
ax.set_ylabel('Solve Time (s)', fontsize=14)
ax.set_xlabel('Number of Similar Paths', fontsize=14)
ax.legend(facecolor='white', fontsize=11)
ax.set_title('$R^2= %.2f$' % r2, fontsize=18)
                         


#%% Let's try this with a log transformation on the times first?
X=pd.DataFrame(data_array[:,1:4], columns=['Path_Length','Nodes_in_Path', 'Similar_Paths'])
y=pd.DataFrame(np.log(data_array[:,0]), columns=['Solve_Time'])
                                          
regr = linear_model.LinearRegression()
model=regr.fit(X, y)
response = model.predict(X)
r2=model.score(X,y)

print('Log Transform to the y axis first')
print(regr.coef_)  
print(r2) 


plt.style.use('default')
plt.style.use('ggplot')

fig, ax = plt.subplots(figsize=(8, 4))

ax.plot(X_path, response, color='k', label='Regression model')
ax.scatter(X_path, y, edgecolor='k', facecolor='grey', alpha=0.7, label='Sample data')
ax.set_ylabel('Log Solve Time (s)', fontsize=14)
ax.set_xlabel('Path Length', fontsize=14)
ax.legend(facecolor='white', fontsize=11)
ax.set_title('$R^2= %.2f$' % r2, fontsize=18)

fig.tight_layout()

#Similar Paths
X_sim=pd.DataFrame(data_array[:,3], columns=['Similar_Paths'])

regr_sim=linear_model.LinearRegression()
model=regr_path.fit(X_sim,y)
response= model.predict(X_sim)
r2=model.score(X_sim,y)   
print('Number of Similar Paths vs Solve Time')
print(regr_path.coef_)
print(r2)

fig, ax = plt.subplots(figsize=(8, 4))

ax.plot(X_sim, response, color='k', label='Regression model')
ax.scatter(X_sim, y, edgecolor='k', facecolor='grey', alpha=0.7, label='Sample data')
ax.set_ylabel('Log Solve Time (s)', fontsize=14)
ax.set_xlabel('Number of Similar Paths', fontsize=14)
ax.legend(facecolor='white', fontsize=11)
ax.set_title('$R^2= %.2f$' % r2, fontsize=18)

#%% Classification

from sklearn import svm
X=np.transpose(np.vstack([data_array[:,1],data_array[:,3]]))

yclass=np.zeros(shape=(2701,1))
for k in range(0,2701):
    t=data_array[k,0]
    if t>=1:
        yclass[k]=1
    else:
        yclass[k]=0

y0=pd.DataFrame(yclass, columns=['Bad_Good'])
plt.scatter(data_array[:,1],data_array[:,3], c=yclass, cmap=plt.cm.Paired, edgecolors='k')

clf=svm.SVC(kernel='linear', class_weight={1:2})

clf.fit(X,np.ravel(yclass))

ax=plt.gca()
xlim=ax.get_xlim()
ylim=ax.get_ylim()
xx=np.linspace(xlim[0], xlim[1], 30)
yy=np.linspace(ylim[0], ylim[1], 30)
YY, XX =np.meshgrid(yy,xx)
xy=np.vstack([XX.ravel(), YY.ravel()]).T

Z=clf.predict(xy).reshape(XX.shape)

a=ax.contour(XX, YY, Z, colors='k', levels=[0], alpha=0.5, linestyles=['-'])


#%%
#plt.contourf(data_array[:,1], data_array[:,3],data_array[:,0] )
plt.scatter(data_array[:,1], data_array[:,3], c=data_array[:,0], marker=".")
plt.colorbar()
