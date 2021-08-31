# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 11:48:40 2021

@author: spunlag
"""

import networkx as nx
import time
import matplotlib.pyplot as plt
from shortestpath_networkx import create_asymmetric_uncertainty_shortest_path_interdiction_nx
import random
import numpy as np
from pyomo.environ import *
import pickle


#Two dimensions

random.seed(a=631996)

p1=0.25
p2=0.75
q1=random.uniform(0.25*p1,0.75*p1)
q2=random.uniform(0.25*p2,0.75*p2)

p1=0.8
p2=0.8
q1=0.2
q2=0.2

l1=-np.log(p1)
l2=-np.log(p2)
r1=-np.log(q1)+np.log(p1)
r2=-np.log(q2)+np.log(p2)

vdim=2

v1=np.linspace(-1,1,100)
v2=np.linspace(-1,1,100)

#%%
plt.figure(1)
plt.axes()
plt.title('Potential Arc Lengths')
plt.xlabel('Arc 1')
plt.ylabel('Arc 2')

for w in v1:
    for u in v2: 
        if sqrt(w*w+u*u)<=1:
            plt.plot(l1+r1*(1+w),l2+r2*(1+u),'bo')
plt.plot(l1,l2,'r*')
plt.text(l1,l2,'(l1,l2)')
plt.plot(l1+r1,l2+r2,'r*')
plt.text(l1+r1, l2+r2, '(l1+r1, l2+r2)')

            
plt.figure(2)
plt.axes()
plt.title('Potential Arc Probabilities')  
plt.xlabel('Arc 1')
plt.ylabel('Arc 2') 
    
for w in v1:
    for u in v2:
        if (sqrt(w*w + u*u))<=1:
            plt.plot(exp(-(l1+r1*(1+w))),exp(-(l2+r2*(1+u))),'ro')  
            
plt.plot(q1,q2,'b*') 
plt.text(q1,q2,'(q1,q2)')  
plt.plot(p1,p2,'b*')
plt.text(p1,p2,'(p1,p2)')      

#%%
rtwid1= p1-q1
rtwid2= p2-q2
plt.figure(3)
for w in v1:
    for u in v2:
        if (sqrt(w*w + u*u)<=1):
            plt.plot(q1+ rtwid1*(u),q2+rtwid2*(w), 'ro')

plt.plot(q1,q2,'b*')
plt.text(q1,q2,'(q1,q2)')  
plt.plot(p1,p2,'b*')
plt.text(p1,p2,'(p1,p2)')
plt.title('Potential Arc Probabilities')
plt.xlabel('Arc 1')
plt.ylabel('Arc 2')


#%% 
r1=-np.log(q1)+np.log(p1)
r2=-np.log(q2)+np.log(p2)
R11=r1
R22=r2
plt.figure(4)
plt.title('Potential Arc Lengths')
plt.xlabel('Arc 1')
plt.ylabel('Arc 2')
for w in v1:
    for u in v2:
        if (sqrt(w*w + u*u)<=1):
            #r1=(-np.log(q1+(p1-q1)*u)+np.log(p1))/(1+u)
            #r2=(-np.log(q2+(p2-q2)*w)+np.log(p2))/(1+w)
            plt.plot(l1+r1 + R11*w,l2+r2 + R22*u,'bo')
            
plt.plot(l1,l2,'r*')
plt.text(l1,l2,'(l1,l2)')
plt.plot(l1+r1,l2+r2,'r*')
plt.text(l1+r1, l2+r2, '(l1+r1, l2+r2)')  


#%%
plt.figure(5) 
plt.title('Potential Arc Probabilities')
plt.xlabel('Arc 1')
plt.ylabel('Arc 2')
for w in v1:
    for u in v2:
        if (sqrt(w*w + u*u)<=1):
            #r1=(-np.log(q1+(p1-q1)*u)+np.log(p1))/(1+u)
            #r2=(-np.log(q2+(p2-q2)*w)+np.log(p2))/(1+w)
            plt.plot(exp(-(l1+r1 + R11*w)),exp(-(l2+r2 + R22*u)),'ro')  
            
plt.plot(q1,q2,'b*')
plt.text(q1,q2,'(q1,q2)')  
plt.plot(p1,p2,'b*')
plt.text(p1,p2,'(p1,p2)')