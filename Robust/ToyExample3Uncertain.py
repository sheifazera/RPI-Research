# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 08:22:18 2019

@author: She'ifa
"""

# -*- coding: utf-8 -*-
#"""
#Created on Thu Jun 13 10:05:33 2019

#"""

import numpy as np
tol=0.1
mU = 2
nL= 2 

nR = 1

nZ = 1

mR = 1

mZ = 1
 

Px_array=np.identity(mR).reshape(mR,mR)
Px_array=tol*Px_array
Px={}
for i in range(0,mR):
    for j in range(0,mR):
        Px[(i+1),(j+1)]=Px_array[i][j]

Py_array=np.identity(mZ).reshape(mZ,mZ)
Py_array=tol*Py_array
Py={}
for i in range(0,mZ):
    for j in range(0,mZ):
        Py[(i+1),(j+1)]=Py_array[i][j]




AR_array = np.array([0, 6]).reshape(2,1).tolist()
AR={}
for i in range(0,2):
    for j in range(0,1):
        AR[(i+1,j+1)]=AR_array[i][j]
        
        
AZ_array = np.array([7,9]).reshape(2,1).tolist()
AZ={}
for i in range(0,2):
    for j in range(0,1):
        AZ[(i+1,j+1)]=AZ_array[i][j]

BR_array = np.array([5, 10]).reshape(2,1).tolist()
BR={}
for i in range(0,2):
    for j in range(0,1):
        BR[(i+1,j+1)]=BR_array[i][j]

BZ_array = np.array([7,2]).reshape(2,1).tolist()
BZ={}
for i in range(0,2):
    for j in range(0,1):
        BZ[(i+1,j+1)]=BZ_array[i][j]

r_array = np.array([62, 117]).reshape(2,1).tolist()
r={}
for i in range(0,2):
    r[i+1]=r_array[i][0]    
    
cR_array = np.array([20]).reshape(1,1).tolist()
cR={}
for i in range(0,1):
    cR[i+1]=cR_array[i][0]
 
cZ_array =  np.array([-38]).reshape(1,1).tolist()
cZ={}
for i in range(0,1):
    cZ[i+1]=cZ_array[i][0]
dR_array = np.array([1]).reshape(1,1).tolist()
dR={}
for i in range(0,1):
    dR[i+1]=dR_array[i][0]
    
dZ_array =  np.array([42]).reshape(1,1).tolist()
dZ={}
for i in range(0,1):
    dZ[i+1]=dZ_array[i][0]

QR_array = np.array([8,9]).reshape(nL,mR).tolist()
QR={}
for i in range(0,2):
    for j in range(0,1):
        QR[(i+1,j+1)]=QR_array[i][j]
        
QZ_array = np.array([0,0]).reshape(nL,mZ).tolist() #needs to be vertical
QZ={}
for i in range(0,2):
    for j in range(0,1):
        QZ[(i+1,j+1)]=QZ_array[i][j]

PR_array = np.array([2,2]).reshape(nL,nR).tolist()
PR={}
for i in range(0,2):
    for j in range(0,1):
        PR[(i+1,j+1)]=PR_array[i][j]
PZ_array =np.array([8,1]).reshape(nL,nZ).tolist() #needs to be vertical
PZ={}
for i in range(0,2):
    for j in range(0,1):
        PZ[(i+1,j+1)]=PZ_array[i][j]
        
s_array = np.array([53, 28]).reshape(2,1).tolist() #needs to be vertical
s={}
for i in range(0,2):
    s[i+1]=s_array[i][0]
    
    
wR_array = np.array([39]).reshape(1,1).tolist()
wR={}
for i in range(0,1):
    wR[i+1]=wR_array[i][0]
    
wZ_array = np.array([27]).reshape(1,1).tolist()
wZ={}
for i in range(0,1):
    wZ[i+1]=wZ_array[i][0]