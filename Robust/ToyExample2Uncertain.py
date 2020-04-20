# Toy Example 2

import numpy as np
tol=0.1
mU = 2
nL= 2 

nR = 0

nZ = 1

mR = 0

mZ = 1
 
AR = {}
zx=nR
zy=nZ

NR_array=np.identity(nR).reshape(nR,nR)
NR_array=tol*NR_array
NR={}
for i in range(0,nR):
    for j in range(0,nR):
        NR[(i+1),(j+1)]=NR_array[i][j]

NZ_array=np.identity(nZ).reshape(nZ,nZ)
NZ_array=tol*NZ_array
NZ={}
for i in range(0,nZ):
    for j in range(0,nZ):
        NZ[(i+1),(j+1)]=NZ_array[i][j]


MR_array=np.identity(mR).reshape(mR,mR)
MR_array=tol*MR_array
MR={}
for i in range(0,mR):
    for j in range(0,mR):
        MR[(i+1),(j+1)]=MR_array[i][j]

MZ_array=np.identity(mZ).reshape(mZ,mZ)
MZ_array=tol*MZ_array
MZ={}
for i in range(0,mZ):
    for j in range(0,mZ):
        MZ[(i+1),(j+1)]=MZ_array[i][j]

AZ = {}
AZ_array=np.array([-2,1]).reshape(mU,mZ)
for i in range(0,mU):
    for j in range(0,mZ):
        AZ[(i+1,j+1)]=AZ_array[i][j]
                

BR = {}

BZ = {}
BZ_array=np.array([3,1]).reshape(mU,nZ)
for i in range(0,mU):
    for j in range(0,nZ):
        BZ[(i+1,j+1)]=BZ_array[i][j]
        
cR = {}
 
cZ_array = np.array([-1]).reshape(mZ,1)

cZ={}
for i in range(0,mZ):
    cZ[i+1]=cZ_array[i][0]

dR = {}

dZ_array = np.array([-2]).reshape(nZ,1)
dZ={}
for i in range(0,nZ):
    dZ[i+1]=dZ_array[i][0]

PR = {}

PZ_array = np.array([1,1]).reshape(nL,nZ) 
PZ={}
for i in range(0,nL):
    for j in range(0,nZ):
        PZ[(i+1,j+1)]=PZ_array[i][j]

s_array  = np.array([-3,30]).reshape(nL,1) 
s={}
for i in range(0,nL):
    s[i+1]=s_array[i][0]

QR = {}

QZ_array = np.array([-3,3]).reshape(nL,mZ) 
QZ={}
for i in range(0,nL):
    for j in range(0,mZ):
        QZ[(i+1,j+1)]=QZ_array[i][j]

wR = {}

wZ_array = np.array([1]).reshape(nZ,1)
wZ={}
for i in range(0,nZ):
    wZ[i+1]=wZ_array[i][0]

r_array  = np.array([12,14]).reshape(mU,1) 
r={}
for i in range(0,mU):
    r[i+1]=r_array[i][0]

