################################################################################
# Copyright 2022 Avalanche Warning Service Tyrol                               #
################################################################################
# This is free software you can redistribute/modify under the terms of the     #
# GNU Lesser General Public License 3 or later: http://www.gnu.org/licenses    #
################################################################################

# -------- Test data for function ( first iteration step from SNPCAMP data) ------- #
#alp = 38
#calcFEM = 0
#E= [ 1172470.92044471, 605854.839145144, 862746.607687970, 892877.495346702, 962059.456524898, 644887.739186792, 12002649.4584851]
#F = 780
#H = [ 0.0470442419903213, 0.188570973338089, 0.0136325860373963, 0.0118989623794615, 0.00480686559700100, 0.00500000000000000, 0.400000000000000 ]
#msswl = np.NAN
#mute = 1
#rho = [ 125.552093802345, 81.3167989970748, 105, 107.300000000000, 112.300000000000, 85.5000000000000, 281.395031395032 ]
#tau_p = 120
# should return:
#msswl,skierloadE,equiHslab,mssANA = [nan, 628935.0600920586, 0.26595362934226907, 638.4301025263883]
# -------- Test data for function ( first iteration step from SNPCAMP data) ------- #

import numpy as np


def get_S_rb15_v2(H,rho,E,tau_p,alp,mute=1,calcFEM=0, PS=None): 
    """Subroutine to calculate the failure initiation (Reuter et al. 2015)
    - includes ANSYS run for max shear stress at the weak layer
    - 170817 calling get_S_rb15([.1,.1,.1,.005,.4],[200,200,200,200,200],[3,3,3,3,3]*1e6,1e5,38)
        and poi=.49 and tuning: F/1.8 and Wski/2 reproduces analytical solution by 99.5%
    - 170817 routine reproduced MAIN_ini_multi_pene.m (watch for WL-thickness!)

    Args:
        layer thickness (H)
        denisty (rho)
        elastic modulus (E)
        strength of WL
        (tau_p)
        slope angle (alp)
        suppress output (mute=1)
        analytic solutions only (calcFEM=0)

       optional: penetration depth (PS)
    Returns:
        msswl
        skierloadE
        equiHslab
        mssANA 
    """

    if PS == None:
        pene = 0
    else:
        pene = 1
    
    """Constants"""
    # poi  = 0.25
    Wski = 0.2
    F    = 780

    """OPTION 1: compact the snow in the slab within the penetration depth"""
    if pene:
        print('not implemented... line33...52')

    """OPTION 2: /wo considering penetration depth"""
    ###otherwise we use H, rho, E as are
    
    if calcFEM:
        print('FEM not implemented', 'matlabcode get_S_rb15_v2 line 54...85')
        ###MODEL stress at the depth of the weak layer
    else:
        msswl = np.NaN
        if not mute:
            print('--- Entering failure ini analyt. sol. --- ')
    
    """Analytic solution"""
    ### Account for layering /w equivalent slab thickness (Monti et al. 2016)   
    skierloadE = (np.sum(np.array(H[0:-2])*np.array(E[0:-2])**0.33)/ np.sum(H[0:-2]))**3
    equiHslab= np.sum(H[0:-2]) #* (skierloadE/E(end-1))^.33;
    
    ### Analytic (surface load) McClung & Sz 99 (see maxshstr_anal_fem.m)
    lineload = F/ 1.8   #force divided by length of skis (def=780/1.8)
    B = lineload / (Wski/2) / (2*np.pi)* np.sin(np.deg2rad(alp))
    A = - lineload / (Wski/2) / (2*np.pi)* np.cos(np.deg2rad(alp))
    theta_deg = np.linspace(0.1,89,101)
    theta = np.deg2rad(theta_deg)
    
    deltatauxy = B * ( np.sin(theta) * np.cos(theta) - np.sin( np.arctan(np.tan(theta)+Wski*2/equiHslab)) \
                          * np.cos( np.arctan(np.tan(theta)+Wski*2/equiHslab)) \
                      ) \
                 - B * ( theta - np.arctan(np.tan(theta) + Wski*2 / equiHslab))   \
                 - 0.5 * A* (  (np.cos(theta)**2 - np.sin(theta)**2) \
                              - ( np.cos(np.arctan( np.tan(theta) + Wski*2/equiHslab))**2 \
                                  - np.sin( np.arctan(np.tan(theta)+Wski*2/ equiHslab))**2 \
                                ) \
                            )
    mssANA = np.max(deltatauxy)
    
    """Gaume Reuter 17 Skier crack length"""
    # F=500;  %default 780 [N][[
    #### idea: get F from: msswl==Boussineq solution (2nd term in next line)
    # 
    # myfun = @(x) tau_g +   2*F/(hslab*cosd(alp)*3.14)*(sind(x))^2*sind(x+alp)*cosd(x)   - tau_p;
    # S1=fzero(myfun,20); %search myfun(x)=0 for values of x starting with initial value 20°
    # S2=fzero(myfun,70);  %search myfun(x)=0 for values of x starting with initial value 70°
    # skiercracklength=hslab*cosd(alp)*(1/tand(S1)-1/tand(S2));

    if not mute:
        print('      Analytic solution  |  FE solution')
        print('max. stress at WL = ', str(round(mssANA)),'  |  ', str(round(msswl)))
        print('failure init crit = ', str(round(tau_p/mssANA*10)/10),'  |  ', str(round(tau_p/msswl*10)/10))

    return msswl,skierloadE,equiHslab,mssANA