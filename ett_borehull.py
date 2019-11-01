#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 13:02:26 2019

@author: joel

Temperaturprofilen i borehullveggen for et felt bestående av ett borehull er
evauert ved å bruke verktøyet pygfunction. 

Denne koden returenrer
    1) Temperaturprofilen i borehullveggen for ett borefelt bestående av ett
       borehull. 

"""

from __future__ import absolute_import, division, print_function

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator
import pygfunction as gt

def main():
    
    #------------------------------------------------------------------------------
    # Parametere i sumuleringen
    #------------------------------------------------------------------------------
    
    # Antall plottpunkter
    n = 100
    
    # Dimensjonering av borehull
    H = 160         # Borehulets lengde (m)
    D = 5           # Borehullets avstand fra overflaten (m)
    r_b = 0.05      # Borehullets radius (m)
    x = 0           # Borehullhodets posisjon langs x - aksen
    y = 0           # Borehullhodets posisjon langs y - aksen
    
    # Termiske egenskaper
    T_g = 6                         # uforstyrret grunntemperatur['C]
    k_g = 1                         # konduktivitet i grunnen [W/mK]
    a_g = 10**-6                    # diffusivitet i grunnen [m^2/s]
    q = 2*np.pi*k_g                 # uniform varmeuttømming fra grunnen [W/m]
    
    # Geometrisk ekspanderende tidsvektor
    dt = 3600                       # Tidssteg [sec]
    aar = 60*60*24*365              # 1. år [sec]
    tmaks = aar*3000                 # Simuleringsperiode [sec]
    Nt = n                          # Ant. tidssteg
    tid = gt.utilities.time_geometric(dt,tmaks,Nt)
    t_s = H**2/(9*a_g)              # Borefeltets karakteristiske tid
    tid_skala = np.log(tid/t_s)     # dimensjonsløs logaritmisk tidsskala
    
    #------------------------------------------------------------------------------
    # Definerer borehullet
    #------------------------------------------------------------------------------
    
    b =[gt.boreholes.Borehole(H, D, r_b, x, y)]
    
    # Visualisering av borehull
    #gt.boreholes.visualize_field(b)
    
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
    ax1.set_xlim([-10.0, 5.0])
    ax1.set_ylim([0., 20.])
    # ticks
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    #grid
    ax1.grid(color = 'k', linestyle='--', linewidth=0.1)
    # Justering av plottvindu
    plt.tight_layout()
    #------------------------------------------------------------------------------
    # Evaluerer g-funksjonen for uniform varmeuttømming fra grunnen
    #------------------------------------------------------------------------------
    for field in [b]:
        gfunk = gt.gfunction.uniform_heat_extraction(b,tid, a_g, disp=True)
        
        # Temperatur i borehullveggen
        T_b1 = T_g + q/(2*np.pi*k_g)*gfunk
        
        #Tegner temperaturprofilen
        ax1.plot(tid_skala, T_b1, 'k-', lw=1.5)
        plt.show()
    
    #------------------------------------------------------------------------------
    # Graf som viser hvilken dag verdien av log(tid/t_s) svarer til
    #------------------------------------------------------------------------------
    
    plt.rc('figure')
    fig = plt.figure()
    ax3 = fig.add_subplot(111)
    # Navn på aksene
    ax3.set_xlabel(r'$ln(t/t_s)$')
    ax3.set_ylabel(r'Dager')
    # aksegrenser
    ax3.set_xlim([-10.0, -4.])
    ax3.set_ylim([0., 350.])
    # ticks
    ax3.xaxis.set_minor_locator(AutoMinorLocator())
    ax3.yaxis.set_minor_locator(AutoMinorLocator())
    #grid
    ax3.grid(color = 'k', linestyle='--', linewidth=0.1)
    # Justering av plottvindu
    plt.tight_layout()
    
    #Tegner graf
    c = np.arange(-10,-4,step=0.2)
    w = np.zeros(len(c))
    for i in range(len(c)):
        w[i] = t_s*np.exp(c[i])
        w[i] = w[i]/(60*60*24)
    
    ax3.plot(c,w)
    ax3.grid(color = 'k', linestyle='--', linewidth=0.1)
    plt.show()
# Main function
if __name__ == '__main__':
    main()
    
    






