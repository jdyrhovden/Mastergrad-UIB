#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 13:50:55 2019

@author: joel
En sammenligning av pygfunction, den uendelige linjesluk-løsningen og den
den uendelige sylindersluk-løsningen for et borefelt bestående av ett borehull 
er gjort. Det er antatt uniform og kontinuerlig varmeuttak langs hele borehullet.

Den uendelige linjeslukløsningen er hentet fra Carslaw og Jaeger sin løsning
fra 1947. 

Koden retunerer: 
    1) Temperaturprofilen ved borehullveggen for de tre modellene beskrevet 
       ovenfor, for et borefelt bestående av ett borehull.
    2) Forskjellen i temperatur ved borehullveggen for de ulike modellene.
    3) En graf som viser hvilket år den dimensjonsløse logaritmiske tidsskalaen
       viser til. 
"""

from __future__ import absolute_import, division, print_function

from matplotlib.legend_handler import HandlerLine2D
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator
import pygfunction as gt
from scipy import integrate
import math as m
import scipy.special as sp
import timeit

# tidtaker 
start = timeit.default_timer()
 
#------------------------------------------------------------------------------
# Parametere i simuleringen
#------------------------------------------------------------------------------

# Antall plottpunkter
n = 100

# Dimensjonering av borehull
H = 160         # Borehulets lengde (m)
D = 5           # Borehullets avstand fra overflaten (m)
r_b = 0.1       # Borehullets radius (m)
x = 0           # Borehullhodets posisjon langs x - aksen
y = 0           # Borehullhodets posisjon langs y - aksen

#

# egenskaper til grunnen
T_g = 6                         # uforstyrret grunntemperatur['C]
k_g = 1                         # konduktivitet i grunnen [W/mK]
a_g = 10**-6                    # diffusivitet i grunnen [m^2/s]
q = 2*np.pi*k_g                 # uniform varmeuttømming fra grunnen [W/m]
r = 20                          # radiell størrelse til grunnen [m]

# Geometrisk ekspanderende tidsvektor
dt = 3600                       # Tidssteg [sec]
aar = 60*60*24*365              # 1. år [sec]
tmaks = aar*3000                # Simuleringsperiode [sec]
Nt = n                          # Ant. tidssteg
tid = gt.utilities.time_geometric(dt,tmaks,Nt)
t_s = H**2/(9*a_g)              # Borefeltets karakteristiske tid
tid_skala = np.log(tid/t_s)     # dimensjonsløs logaritmisk tidsskala

# Tabell med simuleringsparametere
print()
print('Følgende parametere er brukt i simuleringen:')
print()
from tabulate import tabulate
print(tabulate([['q', r'$2\pi \lambda$'], ['H', '160']], 
               headers=['Parameter', 'Verdi'], tablefmt='orgtbl'))
print()

#------------------------------------------------------------------------------
# Definerer borehullet
#------------------------------------------------------------------------------

b =[gt.boreholes.Borehole(H, D, r_b, x, y)]

try:
    # Visualisering av borehull
    gt.boreholes.visualize_field(b)
except Exception:
    pass

#------------------------------------------------------------------------------
# Definerer koordinatsystem
#------------------------------------------------------------------------------

plt.rc('figure')
fig = plt.figure()
ax1 = fig.add_subplot(111)
# Navn på aksene
ax1.set_xlabel(r'$ln(t/t_s)$')
ax1.set_ylabel(r'Temp. ved borehullvegg[($^\circ$C)]')
# aksegrenser
ax1.set_xlim([-14.0, 5])
ax1.set_ylim([0.,30.])
# ticks
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())
#grid
ax1.grid(color = 'k', linestyle='--', linewidth=0.1)
# Justering av plottvindu
plt.tight_layout()

#------------------------------------------------------------------------------
# Beregner g-funksjonen og temperaturprofilen ved borehullveggen
#------------------------------------------------------------------------------

for field in [b]:
    gfunk = gt.gfunction.uniform_heat_extraction(b,tid, a_g, disp=True)
    
    # Temperatur i borehullveggen
    T_b_pyg = T_g + q/(2*np.pi*k_g)*gfunk         

#------------------------------------------------------------------------------
# Evaluerer den uendelige linjesluk - løsningen
#------------------------------------------------------------------------------
   
# Diverse matriser
s = (n,2)
u = np.zeros(n)             # Nedre grense i integralet i ULS-løsningen
I_int = np.zeros(s)         # Verdier av integralet i ULS-løsningen
T_b_uls = np.zeros(s)       # Temperatur ved borehullveggen 

# Evaluerer ULS-løsnningen
for i in range(0,n):  
    # verdier for u=r_b^2/4a_gt
    u[i] = r_b**2/(4*a_g*tid[i])
    
    # integrasjon av integralet i ILS-løsningen
    I = lambda x: (np.exp(-x))/x
    I_int[i] = integrate.quad(I,u[i],m.inf)
    
    #verdier for temperaturen ved borehullveggen
    T_b_uls[i] = T_g + (q/(2*m.pi*k_g))*I_int[i] 

#------------------------------------------------------------------------------
# Evaluerer sylindersluk-løsningen
#------------------------------------------------------------------------------

Fo = np.zeros(n)        # Fourier-tallet
G = np.zeros(s)         # Sylindersluk-funksjonen
T_b_uss = np.zeros(s)   # Temperatur ved borehullveggen for uss løsningen 

for i in range(n):
    # Fouriertallet
    Fo[i] = (a_g/(r_b**2))*tid[i]
    
    # Beregner integralet i sylindersluk-funksjonen
    G[i] = integrate.quad(lambda h:((np.exp(-(h**2)*Fo[i]))-1)*((sp.j0(h)*sp.y1(h) 
    - sp.y0(h)*sp.j1(h))/((h**2)*((sp.j1(h)**2)+(sp.y1(h)**2)))),0,m.inf)
    
    # Temperatur ved borehullveggen
    T_b_uss[i] = T_g + q/(k_g*(np.pi**2))*G[i]
    


#------------------------------------------------------------------------------
# Temperaturprofilen ved borehullveggen for de ulike modellene
#------------------------------------------------------------------------------

line1, = ax1.plot(tid_skala, T_b_pyg,'k-', lw=1.5, label='Pygfunction' )
line2, = ax1.plot(tid_skala, T_b_uls[:,0], 'r-', lw=1.5, label='u_linjesluk')
line3, = ax1.plot(tid_skala, T_b_uss[:,0], 'y-', lw=1.5, label='u_sylindersluk')

plt.legend(handler_map={line3: HandlerLine2D(numpoints=4)}, loc=4)
plt.show()

#------------------------------------------------------------------------------
# Temperaturdifferanse ved borehullveggen for de ulike modellene
#------------------------------------------------------------------------------

plt.rc('figure')
fig = plt.figure()
ax2 = fig.add_subplot(111)
# Tittel
ax2.set_title('Temperaturdifferanse ved borehullveggen')
# Navn på aksene
ax2.set_xlabel(r'$ln(t/t_s)$')
ax2.set_ylabel(r'Temperatur [($^\circ$C)]')
# aksegrenser
ax2.set_xlim([-14.0, 5])
ax2.set_ylim([0.,1.])
# ticks
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())
#grid
ax2.grid(color = 'k', linestyle='--', linewidth=0.1)
# Justering av plottvindu
plt.tight_layout()

# pygfunction og sylindersluk-modellen
D_pyg_uss = np.zeros(n)
D_pyg_uls = np.zeros(n)
T_b_uls = T_b_uls[:,0]
T_b_uss = T_b_uss[:,0]

for i in range(n):
    D_pyg_uss[i] = -(T_b_pyg[i] - T_b_uss[i])
    D_pyg_uls[i] = T_b_uls[i] - T_b_pyg[i]

line1, = ax2.plot(tid_skala,D_pyg_uss, 'y-', lw=1.5, label='uss-pyg')
line2, = ax2.plot(tid_skala,D_pyg_uls, 'b-', lw=1.5, label='uls-pyg')


plt.legend(handler_map={line2: HandlerLine2D(numpoints=4)}, loc=1)
plt.show()

        

#------------------------------------------------------------------------------
# Graf som viser hvilket år verdien av log(tid/t_s) svarer til
#------------------------------------------------------------------------------

plt.rc('figure')
fig = plt.figure()
ax3 = fig.add_subplot(111)
# Navn på aksene
ax3.set_xlabel(r'$ln(t/t_s)$')
ax3.set_ylabel(r'År')
# aksegrenser
ax3.set_xlim([-7.0, -2.])
ax3.set_ylim([0., 10])
# ticks
ax3.xaxis.set_minor_locator(AutoMinorLocator())
ax3.yaxis.set_minor_locator(AutoMinorLocator())
#grid
ax3.grid(color = 'k', linestyle='--', linewidth=0.1)
# Justering av plottvindu
plt.tight_layout()

#Tegner graf
c = np.arange(-10, 3,step=0.2)
w = np.zeros(len(c))
for i in range(len(c)):
    w[i] = t_s*np.exp(c[i])
    w[i] = w[i]/(60*60*24*365)

ax3.plot(c,w)
ax3.grid(color = 'k', linestyle='--', linewidth=0.1)
plt.show()


# tidtaker
stop = timeit.default_timer()
print('Total beregningstid: ', stop - start, 'sec')
    

    
    








