# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 09:44:08 2021

@author: spunlag
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle

plt.figure()
plt.axes()
ax = plt.gca()

ellipse = Ellipse(xy=(0.5, 0.4), width=1, height=0.8, 
                        edgecolor='r', fc='r', lw=2)

rect = Rectangle((0, 0),  1, 0.8, linewidth=1, edgecolor='k', facecolor='None')
ax.add_patch(ellipse)
ax.add_patch(rect)
ax.plot(0.5,0.4, 'ko')
plt.text(0.3,0.43,r'($\ell_1 +\bar{r_1}$,$\ell_2 + \bar{r_2}$ )=(0.5,0.4)')  
ax.plot(1,0.4,'ko')
plt.text(0.67, 0.35, r'($\ell_1$,$\ell_2 + \bar{r_2}$) = (1,0.4)')
ax.plot(0.5, 0.8, 'ko')
plt.text(0.35, 0.73, r'($\ell_1 +\bar{r_1}$,$\ell_2$ )=(0.5,0.8)')
