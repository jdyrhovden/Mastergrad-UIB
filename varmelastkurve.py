#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 19:31:42 2019

@author: joel

To funksjoner defineres i dette skriptet: 
    1) Heavysidefunksjonen
    2) Varmelastkurve for stykkevis konstante varmelaster

"""
from __future__ import absolute_import, division, print_function

import matplotlib.pyplot as plt
import numpy as np


#definerer hevysidefunksjonen
def hf(x):
    
    """
    .....
    
    Definisjonen av Heavysidefunksjonen
    
    Parameter
    --------
        x : inputvariabel i funksjonen
    
    Returnerer
    ---------
        Verdien til Heavysidefunksjonen
    """
    if x == 0:
        return 1
    return 0 if x < 0 else 1



def varmelastkurve(per: int, per_enhet: str, q: list): 
    """
    .....
    
    Tegner varmelastkurven for stykkevis konstante varmelaster. 
    
    Paramtetere
    -----------
    
        per       : Perioden som evalueres (måneder).
        per_enhet : Enheten til perioden som evalueres.
        q         : liste som inneholder tall som viser hvordan 
                    varmebelastningen øker/minker per tidssteg. 
    
    Returnerer
    ----------
        
        Grafen til den varierende varmelasten q(t).
    """
    
    #tidssteg
    t_steg = len(q)+1
    t = np.linspace(0,per,t_steg)
    
    #enhet x - aksen
    en = per_enhet
    
    #plotpunkt
    n = 1000
    x = np.linspace(0,per,n)
    
    #Evaluerer q(t)
    q_sum = np.full(n,q[0]) #summen av varmebelastninger
    q_sum[0] = q[0]
    for j in range (1,len(q)):
        for i in range(0,len(x)):
            q_sum[i] = q_sum[i] + q[j]*hf(x[i]-t[j])
    
    #Definerer koordinatsystem
    plt.xlabel("Evalueringsperiode {}".format(en))
    plt.ylabel('Varmelast [W/m]')
    plt.title('Dekomponering av variabel varmelast q(t)')
    plt.grid(color = 'k', linestyle='--', linewidth=0.1)
    
    #Tegner de individuelle varmelastene
    for i in range(len(q)):
        plt.hlines(q[i],xmin=t[i], xmax=t[-1], color='r', linestyle='--', 
                   linewidth = 0.5)
        plt.plot(t[i],q[i], 'o')
        
    #returnerer grafen til den totale varmelasten i hele evalueringsperioden                
    return [plt.plot(x,q_sum, 'k'), plt.show(), 
            print("\n Følgende tidssteg evalueres ({}):".format(en), t, '\n')] 
    
            
            
            


    
   







    








    


