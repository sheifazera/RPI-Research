# -*- coding: utf-8 -*-
'''
@author: She'ifa
'''


import numpy as np
delta=0.01
nu=2
mu=2
ml=2
nl=2
G={}
H={}
q={}
A={}
B={}
b={}
cx={}
cy={}
d={}

cx_array = np.array([-1,-1]).reshape(nu,1)

for i in range(0,nu):
    cx[i+1]=cx_array[i][0]

cy_array = np.array([-1,-1]).reshape(nl,1)

for i in range(0,nl):
    cy[i+1]=cy_array[i][0]
    
    
A_array=np.array([[1,2],[3,4]]).reshape(ml,nu)  
for i in range(0,nu):
    for j in range(0,ml):
        A[(i+1,j+1)]=A_array[i][j]
        
B_array=np.array([[1,2],[3,4]]).reshape(ml,nl)  
for i in range(0,nl):
    for j in range(0,ml):
        B[(i+1,j+1)]=B_array[i][j]
                
        
b_array=np.array([0,0]).reshape(ml,1)

for i in range(0,ml):
    b[i+1]=b_array[i][0]
   
    
H_array=np.array([[1,2],[3,4]]).reshape(mu,nl)  
for i in range(0,mu):
    for j in range(0,nl):
        H[(i+1,j+1)]=H_array[i][j]
        
G_array=np.array([[1,2],[3,4]]).reshape(mu,nu)  
for i in range(0,nu):
    for j in range(0,mu):
        G[(i+1,j+1)]=G_array[i][j]
                
        
q_array=np.array([0,0]).reshape(mu,1)

for i in range(0,mu):
    q[i+1]=q_array[i][0]
        
d_array=np.array([1,1]).reshape(nl,1)

for i in range(0,nl):
    d[i+1]=d_array[i][0]