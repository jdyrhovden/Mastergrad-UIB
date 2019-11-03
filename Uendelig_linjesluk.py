#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 09:06:01 2019

@author: joel

Den uendelige linjesluk - løsningen gitt av Carslaw og Jaeger (1947) er 
evaluert for en konstant varmeuttømming fra grunnnen, og tar for seg 
to tilfeller:
    1) Tiden holdes konstant og radius fra linjesluket varierer.
    2) Radius holdes konstant (=radius til borehullveggen) og tiden varierer.

Denne koden returnerer
    1) Temperaturprofilen i grunnen i en avstand r[m] fra linjesluket for
       perioder på: 1 dag, 1 måned, 6 månder og 1 år. 
    2) Temperaturprofilen ved borehullveggen for en periode på 1 år. 
"""

from scipy import integrate
from matplotlib.ticker import AutoMinorLocator
from matplotlib.legend_handler import HandlerLine2D
import matplotlib.pyplot as plt
import numpy as np
import math as m
import timeit

def main():
    start = timeit.default_timer()
    
    #------------------------------------------------------------------------------
    # Parametere i sumuleringen
    #------------------------------------------------------------------------------
    
    #Antall plott-punkter
    n = 100
    
    # tidsperioder
    t_1 = 60*60*24                 # 1 dag
    t_2 = t_1*31                   # 1 måned
    t_3 = t_2*6                    # 6 måneder
    t_4 = t_1*365                  # 1 år
    t_5 = t_4*10                   # 10 år
    t = [t_1, t_2, t_3, t_4, t_5]  #tidsvektor
    
    # (Termiske) egenskaper 
    T_g = 6                        # uforstyrret grunntemperatur['C]
    k_g = 1                        # konduktivitet i grunnen [W/mK]
    a_g = 10**-6                   # diffusivitet i grunnen [m^2/s]
    q = -2*np.pi*k_g               # uniform varmeuttømming fra grunnen [W/m]
    
    # Egenskaper til borehull
    r_b = 0.05                     # Avstand til borehullveggen[m]
      
    # Diverse matriser
    u = np.zeros(n)     # verdier for nedre grense i intergralet i ILS-løsningen
    s = (n,2)
    I_int = np.zeros(s) # Verdier av integralet i ILS-løsningen
    T = np.zeros(s)     # Temperatur i avstand r[m] fra borehullet
    z = (n,len(t))
    T_t = np.zeros(z)   # Temperatur r[m] fra borehullet for ulike tidsperioder
                        
    
    # Verdier som skal evalueres i radiell avstand fra linjesluk
    r_g = 20                        # Avstand fra linjesluket som evalueres [m]
    r = np.linspace(0,r_g,num=n)
    
    #------------------------------------------------------------------------------
    # Definerer koordinatsystem
    #------------------------------------------------------------------------------
    
    plt.rc('figure')
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    # Tittel
    ax1.set_title('Temperaturprofil ved varmeuttømming')
    # Aksenavn
    ax1.set_xlabel(r'Radiell avstand fra linjesluk[m]')
    ax1.set_ylabel(r'Temperatur[($^\circ$C)]')
    # Aksegrenser
    ax1.set_xlim([0., 20.])
    ax1.set_ylim([-2., 7.])
    # ticks
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    # grid
    ax1.grid(color = 'k', linestyle='--', linewidth=0.1)
    # Justerng av plottvindu
    plt.tight_layout()
    
    #------------------------------------------------------------------------------
    # Evaluering av ILS-løsningen, der t=konstant og radius fra borehull varierer
    #------------------------------------------------------------------------------
    
    for j in range(0,len(t)-1):
        
        for i in range(0,n):  
            # verdier for u=r^2/4at
            u[i] = r[i]**2/(4*a_g*t[j])
            
            # integrasjon av integralet i ILS-løsningen
            I = lambda x: (np.exp(-x))/x
            I_int[i] = integrate.quad(I,u[i],m.inf)
            T[i] = T_g + (q/(2*m.pi*k_g))*I_int[i]     
       
        T_t[:,j] = T[:,0] 
       
    # plott av temperaturprofilene 
    line1, = ax1.plot(r, T_t[:,0],'r-', markersize=0.5, label='Etter 1 dag')      
    line2, = ax1.plot(r, T_t[:,1],'b-', markersize=0.5, label = 'Etter 1 måned')  
    line3, = ax1.plot(r, T_t[:,2],'y-', markersize=0.5, label = 'Etter 6 måneder')  
    line4, = ax1.plot(r, T_t[:,3],'k-', markersize=0.5, label= 'Etter 1 år')
    
    plt.legend(handler_map={line4: HandlerLine2D(numpoints=4)}, loc=4)
    
    plt.show()
    
    #------------------------------------------------------------------------------
    # Evaluering av ILS-løsning, der radius er konstant (=r_b) og tiden varierer
    #------------------------------------------------------------------------------
    
    # Tidssteg
    tid = np.linspace(0,t_4, num = n) 
    
    # Diverse matriser
    u1 = np.zeros(n)        # Nedre grense i integralet i ILS-løsningen
    I_int1 = np.zeros(s)    # Verdier av integralet i ILS-løsningen
    T1 = np.zeros(s)        # Temperatur ved borehullveggen som funksjon av tid
    
    # Evaluerer ILS-løsnningen
    for i in range(0,n):  
        # verdier for u=r_b^2/4a_gt
        u1[i] = r_b**2/(4*a_g*tid[i])
        
        # integrasjon av integralet i ILS-løsningen
        I = lambda x: (np.exp(-x))/x
        I_int1[i] = integrate.quad(I,u1[i],m.inf)
        
        #verdier for temperaturen ved borehullveggen
        T1[i] = T_g - (q/(2*m.pi*k_g))*I_int1[i]  
    
    # Definerer koordinatsystem
    plt.rc('figure')
    fig = plt.figure()
    ax2 = fig.add_subplot(111)
    # Tittel
    ax2.set_title('Temperaturprofil ved borehullveggen')
    # Aksenavn
    ax2.set_xlabel(r'Tid(sec)')
    ax2.set_ylabel(r'Temperatur[($^\circ$C)]')
    # Aksegrenser
    ax2.set_xlim([0.,t_4])
    ax2.set_ylim([5., 17.])
    # ticks
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    # grid
    ax2.grid(color = 'k', linestyle='--', linewidth=0.1)
    # Justerng av plottvindu
    plt.tight_layout()     
    
    # plott av temperaturprofil 
    line1, = ax2.plot(tid, T1[:,0],'r-', markersize=0.5, label='Temp. ved radius: r_b')     
    
    plt.legend(handler_map={line4: HandlerLine2D(numpoints=4)}, loc=4)
    plt.show()
    
    stop = timeit.default_timer()
    
    print('Total beregningstid ', stop - start, 'sec')

# Main function
if __name__ == '__main__':
    main()



