# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 09:35:13 2021

@author: spunlag
"""

import matplotlib.pyplot as plt
import numpy as np

x=np.linspace(0,1)

plt.plot(x,x/10,linestyle='solid', color='k')
plt.plot(x,-x/10,linestyle='solid',color='k')
plt.plot(x,x/10 - 0.025, linestyle='dotted')
plt.arrow(0.5, 0.075, -0.25, 0, head_width=0.0125, color='k')
plt.text(0.25,0.085,'Objective f(x,y)')
plt.plot(0,0, 'co')
plt.plot(0.125, -0.0125, 'ro')
plt.fill_between(x, x/10, -x/10, alpha=0.2)