# -*- coding: utf-8 -*-
#"""
#Created on Wed Jun 12 13:01:02 2019

#@author: spunlag
#"""

import numpy as np

mU = 0
nL= 4 

nR = 0

nZ = 1

mR = 0

mZ = 1
 
AR = {}
AZ = {}
BR = {}
BZ = {}
 
cR = {}
 
cZ_array = np.array([-1]).reshape(mZ,1)

cZ={}
for i in range(0,mZ):
    cZ[i+1]=cZ_array[i][0]

dR = {}

dZ_array = np.array([-10]).reshape(nZ,1)
dZ={}
for i in range(0,nZ):
    dZ[i+1]=dZ_array[i][0]

PR = {}

PZ_array = np.array([20, 2, -1, -10]).reshape(nL,nZ) 
PZ={}
for i in range(0,nL):
    for j in range(0,nZ):
        PZ[(i+1,j+1)]=PZ_array[i][j]

s_array  = np.array([30, 10, 15, -15]).reshape(nL,1) 
s={}
for i in range(0,nL):
    s[i+1]=s_array[i][0]

QR = {}

QZ_array = np.array([-25, 1, 2, -2]).reshape(nL,mZ) 
QZ={}
for i in range(0,nL):
    for j in range(0,mZ):
        QZ[(i+1,j+1)]=QZ_array[i][j]

wR = {}

wZ_array = np.array([-1]).reshape(nZ,1)
wZ={}
for i in range(0,nZ):
    wZ[i+1]=wZ_array[i][0]

r =  {}

