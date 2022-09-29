#%% solver to get critical crack length (Reuter et al. 2015) and analy. sol.
#% i.e. solving the equation presented in Schweizer et al. 2011, which was
#% derived from beam bending eq after Heierli 2008
#% can be used with get_Ebulk.m (for effect. modulus) which includes FE-based calibration after van Herwijnen et al. 2016
#% also incl. analytic solutions of Sigrist 2006 (his diss., before Heierli)
#% also incl. Gaume 2017 (heuristic from DEM results) 
#% mute=1: stops notifications

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import e
import datetime
from fu_tau_p_CoJ15 import fu_tau_p_CoJ15
from dateutil.relativedelta import relativedelta
import get_s_rb15_v2 as pro_rb15
import sympy as sym

# laws and constants:


def Ela(x): # Scapozza
    ela =1.8*10**5 * np.exp(x/67)
    return ela

def Heierlifun ( r,E, H, alp, rho):
    # %Heierli's formula for total energy of a crack 
    gamma = 1
    eta = (12/15)**0.5
    g = 9.81
    tau = - rho * g * H * np.sin(alp/180*np.pi)
    sig = - rho * g * H * np.cos(alp/180*np.pi)

    w = 3/8 * eta**2 * H /E * tau**2 + \
        (0.5 * np.pi * gamma / E + 3/4 * eta /E) * tau**2 * r + \
        3/2 * eta **2 / E * sig * tau * r + \
        0.5 * np.pi * gamma / E * sig**2 * r + \
        0.5 / E /H * tau**2 * r**2 + \
        9/4 * eta / E / H * sig * tau * r**2 + \
        3/2 * eta**2 / E / H * (sig**2) * (r**2) + \
        3/2 * eta / E / (H**2) * (sig**2) * (r**3) + \
        3/2 / E / (H**3)* (sig**2) * (r**4)
    return w

def get_ac_vh16_v2(Eslab,rhoslab,hslab,wf,alp,  mute =1,  rho_wl = -999 ,tau_p= -999 ):
    #% solver to get critical crack length (Reuter et al. 2015) and analy. sol.
    # i.e. solving the equation presented in Schweizer et al. 2011, which was
    # derived from beam bending eq after Heierli 2008
    # can be used with get_Ebulk.m (for effect. modulus) which includes FE-based calibration after van Herwijnen et al. 2016
    #a also incl. analytic solutions of Sigrist 2006 (his diss., before Heierli)
    # also incl. Gaume 2017 (heuristic from DEM results) 
    # mute=1: stops notifications
    
    if not mute :
        print ('--- Entering ac calc /w calibration of vh16 ---')
        
    #%% solve Heierli's eq for critical cut length (r)
    r = sym.Symbol('r', real = True) # just real part of solution
    
    A = sym.solve(Heierlifun(r,Eslab,hslab,alp,np.round(rhoslab))-wf)
    ac = []
    for ele in A:
        if ele > 0: # positive solution
            ac.append(ele)
    if not ac:
        ac.append(np.NaN)
    
    if not mute:
        print('ac = ', str(round(ac*100)/100), 'm;    wf = ', str(round(wf*100)/100),'J/m^2') 
    
    # %% get other analytic solutions too
    tau_g = rhoslab * 9.81 * hslab * np.sin(alp*np.pi/180)
    sigma_g = rhoslab * 9.81 * hslab * np.cos(alp*np.pi/180)
    
    #% Sigrist 2006
    ac_si06 = ( 2 * Eslab * hslab**3 * wf / 3 / (tau_g**2 + sigma_g**2) )**0.25 #  % note that tau^2+sigma^2 = (rhoslab*g*hslab)^2
    
    
    if rho_wl== -999:
        # Gaume 2017 (bending and tension terms)
        ac_ga17 = np.NaN
        if not mute:
            print('data missing to calc ac after Gaume17')
    else:    
        Hwl = 0.01
        Gwl = Ela(rho_wl)/2
        lam = (Eslab / Gwl * hslab * Hwl)**0.5
        ac_ga17 = lam * (-tau_g+ np.emath.sqrt(tau_g**2 + 2*sigma_g*(tau_p-tau_g))) / sigma_g

    
    
    return ac[0], wf, ac_si06, ac_ga17
    
    #return ac, wf, ac_si06, ac_ga17

