#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 09:36:59 2019

@author: joel

Besselfunksjonene av 1. og 2. orden er en del av sylindersluk-modellen for 
varmetransporten mellom borehull og grunnen. 


Denne koden returnerer grafen til Besselfunksjonene av 1. og 2. orden. 

"""

import scipy.special as special
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.legend_handler import HandlerLine2D

#definerer koordinatsystem
plt.rc('figure')
fig = plt.figure()
ax1 = fig.add_subplot(111)
# Navn på aksene
ax1.set_xlabel(r'x')
ax1.set_ylabel(r'Besselfunksjoner')
# Tittel
ax1.set_title('Besselfunksjoner av 1. og 2. orden, J=1.orden, Y=2.orden')
# aksegrenser
ax1.set_xlim([0., 20.])
ax1.set_ylim([-1, 1.1])
# ticks
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())
#grid
ax1.grid(color = 'k', linestyle='--', linewidth=0.1)
# Justering av plottvindu
plt.tight_layout()

# x - verdier
x = np.linspace(0,20,num=500)

# Bessel-integral av 1. orden
J0 = special.j0(x)
J1 = special.j1(x)

# Bessel-integral av 2. orden
Y0 = special.y0(x)
Y1 = special.y1(x)

# Tegner Besselfunksjoner av 1. og 2. orden
line1, = ax1.plot(x,J0, 'r-', markersize=0.5, label='J0')
line2, = ax1.plot(x,J1, 'k-', markersize=0.5, label='J1')
line3, = ax1.plot(x,Y0, 'y-', markersize=0.5, label='Y0')
line4, = ax1.plot(x,Y1, 'b-', markersize=0.5, label='Y1')

plt.legend(handler_map={line4: HandlerLine2D(numpoints=4)}, loc=1)
plt.show()





















