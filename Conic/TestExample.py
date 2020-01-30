# -*- coding: utf-8 -*-
'''
Created on Thu Oct 17 14:40:49 2019

@author: She'ifa
'''


import numpy as np

m=2
n=2
A={}
b={}
c={}
P={}

c_array = np.array([-1,-1]).reshape(n,1)

for i in range(0,n):
    c[i+1]=c_array[i][0]
    
    
A_array=np.array([[1,2],[3,4]]).reshape(m,n)  
for i in range(0,n):
    for j in range(0,m):
        A[(i+1,j+1)]=A_array[i][j]

b_array=np.array([0,0]).reshape(m,1)

for i in range(0,m):
    b[i+1]=b_array[i][0]
    
P_array=np.array([[1,0],[0,1]]).reshape(m,n)  
for i in range(0,n):
    for j in range(0,m):
        P[(i+1,j+1)]=P_array[i][j]