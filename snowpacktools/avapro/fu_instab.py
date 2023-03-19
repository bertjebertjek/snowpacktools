################################################################################
# Copyright 2022 Avalanche Warning Service Tyrol                               #
################################################################################
# This is free software you can redistribute/modify under the terms of the     #
# GNU Lesser General Public License 3 or later: http://www.gnu.org/licenses    #
################################################################################

import numpy as np
import pandas as pd

# from math import e
import datetime
from dateutil.relativedelta import relativedelta

import matplotlib.pyplot as plt

from snowpacktools.avapro import get_s_rb15_v2 as pro_rb15
from snowpacktools.avapro import get_ac_vh16_v2 as pro_vh16
from snowpacktools.avapro import fu_tau_p_CoJ15

# laws and constants:
ini_rc = 0.2
Ltot = 1.2

def Ela(x): # Scapozza
    ela =1.8*10**5 * np.exp(x/67)
    return ela

def TSjj(x): # %Jamieson Johnston
    tsjj = 58300 * (x/917)**2.65
    return tsjj

def SSjj(dens):# strength temporal evolution law, actually strength=cohesion+internalfriction (Gaume etal 14, JGR) 
    ssjj = 18500 * (dens/917)**2.11
    return ssjj

def SScj(dys, c_0): # %(average:^.26, other see ref in Conlan&Jamieson 2015, CGJ)
    sscj = c_0 * dys**0.3
    return sscj

def CERR ( ac, failurestress, Emod): # %(Griffith's eq. or use: Gaume etal 2014, JGR, eq. B7)
    cerr = 1.12**2 *failurestress**2 * ac * np.pi / Emod  
    return cerr

threshprecstab = 2.7 #  %critical value for natural stability (see ConwayWilbur99 eq 6, theort. default: 1)


def aggregate_layers(hh, rho_, trsh, plotit=0, mute=0):
    """ Function to aggregate layers
    - Layers with a difference between them which is smaller than the threshold are averaged and 
    aggregated into one layer, hence reducing the number of layers
    - INPUT:
        h: thickness
        rho: density
        trsh: threshold (def: .1 would be 10%)
    """

    if not mute :
        print ('--- Entering aggregate_layers ---')
        print ('averaging densities up to ', str(trsh*100),' difference')
        
    xx = hh.copy()
    yy = rho_.copy()

    di = (yy[0:-1].values - yy[1:].values) / yy[1:].values
    idn = np.where( abs(di) < trsh)[0]
    # % ct=1; matlab code relict line 21

    while len(idn)> 2 :
        #for kk = 1; matlabrelict line 23
        # %     disp([num2str(ct),'. ','agg: ',num2str(id(1))])
        yy[idn[0]] = (yy[idn[0]] * xx[idn[0]] + yy[idn[0]+1]* xx[idn[0]+1]) / \
                        (xx[idn[0]] + xx[idn[0]+1])
        xx[idn[0]] = xx[idn[0]] + xx[idn[0] + 1]
        xx[idn[0]+1: (len(xx)-1)] = xx[ idn[0]+2 : len(xx) ]
        yy[idn[0]+1: (len(xx)-1)] = yy[ idn[0]+2 : len(xx) ]
        xx = xx[0:-1]
        yy = yy[0:-1]
        di = (yy[0:-1].values - yy[1:].values) / yy[1:].values
        idn = np.where( abs(di) < trsh)[0]
        # %     ct=ct+1; relict from matlabcode line 40
    
    if plotit:
        fig, ax = plt.subplots()
        ax.plot(rho_, -np.cumsum(hh), marker='o',label = 'SNP-output')      
        xb = np.repeat(-np.cumsum(xx), 2)[:-1]
        yb = np.repeat(yy, 2)[1:]
        ax.plot(yb, xb, label ='aggregated layers')
        ax.set_title('A single plot')
        ax.set_xlabel('Density [kg/m^3]')
        ax.set_ylabel('Depth ')
        ax.legend()
        
    if not mute:
        print( str(len(xx),'layers remained'))
    
    return xx, yy

# def sig_dens(hrs):
#     sig_dens = 265 * (1+hrs/24) ** strexpLT 
#     return sig_dens


def fu_precstabindx(burialdate, timewindow, precrat, scmodstep, alp=38, plotit=0 ):
    """strength and natural stability evolution after burial
    - strength fluctuations with natural precipitation loading (see ST and LT)
    - stability after Conway and Wilbur 99 (CRST)
    - calculate precipitation stability index (eq.6): WLshearstrength / load due to precipitation 
    - derive expected time to failure (eq.7):
    ( precstab - threshprecstab )./ (d precstab / d t)
    need to get times with continuous snowfall - this is when strengthening occurs (e.g. estimate from SC model data)
    - all evaluation in hours
    """

    strexpLT = .62 #strength increase exponent for long-term evolution (0.18 from ConlanJamieson15; .63 from Jürg's shear frame data)
    strexpST = .62 #strength increase exponent for short-term, SHOULD DEPEND ON NORMAL LOAD (Szabo and Schneebo, 07)

    # hourly precip shear stress (ConwayWilbur99); 
    precstr = 9.81 * np.cumsum(precrat) * scmodstep * np.cos(alp *np.pi/180)*np.sin(alp*np.pi/180)  #hourly precip shear stress ConwayWilbur99 eq 1
    
    # %hours with continuous snowfall
    #smooth the data
    smoprecstr  = pd.Series(precstr.reset_index(drop=True)).rolling(3, center= True, min_periods=2).mean().to_list()
    #    smoprecstr = precstr.reset_index(drop=True)
    consno = np.zeros(len(smoprecstr))
    iva = smoprecstr[0]
    if np.isnan(iva): # catch nan problem for comparison below
        iva = 0


    ############## INITIALIZATION WITH RTA IS EFFECTED HERE ##############

    for mm in range (1,len(smoprecstr)):
        if smoprecstr[mm] > (iva*1.001): # %watch out, "smooth" makes a tiny increase at the end of the array, so use greater 100.1%
            consno[mm]=consno[mm-1]+1  #% signif. snowfall since last timestep
        else:
            if consno[mm-1] > 0 :   #% no signif. snowfall since last timestep
                consno[mm]=consno[mm-1] - 1/4   #%2: time after snowfall during which strength decreases with dampexp towards longterm trend
        iva = smoprecstr[mm]
    
    consno = consno * scmodstep # %hours with continuous snowfall    
    #% sigLT: long-term strength evolution: SLOW matlab relikt
    #% sig_ini=@(xx)(18.5*1e3*(xx/917).^2.11);   % strength for persistent after JaJo01 
    
    #sig_dens=@(hrs)( 265*(1+hrs/24).^strexpLT );   % strength evolution slow, see  validate_nat_instab_model.m (Jürg's shear frame data from 2015/16)
    Ampl=[]
    for x in timewindow:
        y =  20/100 * np.exp(-x/120) # %20% amplitude dampening in 120h: no more struct changes after many days
        Ampl.append(y)
    dampexp = 1/10   #;  %1/timescale for dampening in 10h: first quick changes then slow changes
    #% WL strength: combine two-speed evolutions
    sigLT = [] # % long-term contribution
    for x in timewindow:
        y = 265 * (1+x/24) ** strexpLT 
        sigLT.append(y)
    

    ###### DEBUG BURIAL DATE FOR INI ######
 
    rrr =  1 + Ampl * (consno) ** strexpST * np.exp(-consno*dampexp)
    sig = sigLT * rrr           #% strength, the sum of LT and ST contribution
    sigST =sigLT * (1-rrr)       #% short-term contribution
    
    ### Stability metrics
    ### precstab = precipitation stability index  ConwayWilbur99 eq 6
    precstab = sig / precstr
    ### Filter out 'inf' values and replace them with np.nan
    precstab = np.where(precstab==np.inf,np.nan,precstab)
    extm2fail = (precstab-1) / np.append(np.NaN, np.abs(np.diff(precstab)))
    extm2fail[extm2fail == 0] = np.NaN

    if plotit:
        datenum = [] # get hourly ticks
        for element in timewindow:
            datetime.timedelta(hours = 1)
            date = burialdate + relativedelta(hours=element)
            datenum.append(date)

        fig, axs = plt.subplots(2,1)
        axs[0].set_xlim(datenum[0], datenum[-1] + datetime.timedelta(hours=19))
        axs[0].plot(datenum, sig, 'r*-', label='sig')
        axs[0].plot(datenum,sigLT,'r--', label='sigLT')
        axs[0].plot(datenum, np.append(0, np.diff(precstr))*10,'--k', label= 'precstr')  
        axs[0].legend()
        axt = axs[0].twinx()
        axt.plot(datenum, np.append(np.NaN, np.diff(consno))>0,'*k', label= 'consno')
        axt.plot(datenum, precstab,'--b' ,label='precstab')
        axt.legend()
        axs[0].set_xlabel('date')
        axs[0].grid(True)
        
        axs[1].set_xlim(datenum[0], datenum[-1] + datetime.timedelta(hours=19))
        axs[1].semilogy(datenum, extm2fail,'.b', label='extm2fail')
        axs[1].semilogy(datenum, extm2fail,'b', label='extm2fail')
        axs[1].set_xlabel('date')
        axs[1].legend()
        axs[1].grid(True)

        fig.tight_layout()
        plt.show()
    
    ### Make sure tm has at least24 entries
    ### Table for sorting: time; precipitation stability index; expected time to failure 
    KKK = np.vstack((timewindow,precstab,extm2fail)).T
    if len(KKK) <=24:# %need at least 24 entries
        KKK = np.concatenate((np.full((24-len(KKK), 3), np.nan), KKK), axis=0) 
    KKK = KKK[-24:,:]
    KKK = KKK[np.argsort(KKK[:, 2])] # %sort by expected time to failure
    
    extm2failMIN24 = np.nanmean(KKK[0][2])  #%expected time to failure, average the 3 minima of the past 24 hours
    precstabMIN24 = np.nanmean(KKK[0][1])   #%min precipitation stability index (corresponding)
    tmcrit = np.min(KKK[0][0])  #%moment when 1st minimum of extm2fail isreached
    tmcrit = np.nanmax(KKK[:,0]) - tmcrit +1 #%hours before drytime when first min of extm2fail was reached

    return precstabMIN24,extm2failMIN24,tmcrit,sig,sigLT,sigST,consno,precstab,extm2fail,precstr


def fu_instab(index,df_prof,df_P,df_met,alp=38,c_0=np.NaN,opt='per', WLopt = 'pap',pre_existing = False ,calcFEM=0,overwrite_scmod_strength=0, winddrift = False, dap_nr = False):
    """Calculate instability of snowprofile
    """
    if WLopt == 'pap':
        WLup =  df_P['papup'][index] # 'papup'
        burialdate = df_P['papburial'][index] #'papburial'
        if pre_existing:
            ini_tau_p = df_P['papINI_tau_p'][index-1]
            c_0 = df_P['papINI_c_0'][index-1]
    elif WLopt == 'nap':
        WLup = df_P['napup'][index] #'napup'
        burialdate = df_P['napburial'][index] # 'napburial'
        if pre_existing:
            ini_tau_p = df_P['napINI_tau_p'][index-1]
            c_0 = df_P['napINI_c_0'][index-1]
    elif WLopt == 'dap':
        WLup = df_P['dapup'][index][dap_nr] #'dap'
        burialdate = df_P['dapburial'][index][dap_nr]
        ini_tau_p = df_P['dapINI_tau_p'][index-1][dap_nr]
        c_0 = df_P['dapINI_c_0'][index-1][dap_nr]
        
    else:
        print('[E]  WL opt not implemented jet -> handling if pap,nap,...')
    if winddrift :
        burialdate=  df_P['dy'][index] #'dy'        

    # read layer file
    # not needed take df_prof
    po = 'height_m' # just for better readability
    rho = 'density'
    shstrength = 'stress in (kPa)'
    h = 'thickness_m'
    lwc = 'lwc'
    gt = 'graintype'
    tcr=np.NaN
    
    if opt == 'windonly' : # % modify profile data, as WL does not exist in snow cover model output
            
        drft = min( [df_P['drft'][index], 1] )  # max thickness of snow drift
        
        slab_prop = {
            'height_m': df_prof[po][0] + drft,
            'density': 200,
            'stress in (kPa)': 1500,
            'thickness_m': drft,
            'graintype': [['RG', 'RG', 'dry']], 
            'lwc': 0
        }
        df_prof = pd.concat([pd.DataFrame(slab_prop), df_prof], ignore_index=True)
        ind = 1
    else:
        ### WL exists in snow cover model output #ind=find(BB.po==WLup);  %WL position
        ind = df_prof.index[df_prof[po]== WLup].values[0]
        
    """Aggregate layers"""
    trsh = .2 # % def: .1, i.e. aggregates if density difference below 10% 
    H, D = aggregate_layers(df_prof[h].iloc[0:(ind)] * np.cos(alp*np.pi/180),df_prof[rho].iloc[0:(ind)], trsh, plotit=0,mute=1) # %last args: plotit, mute
    H_const = pd.Series([0.005,0.4]) 
    H = pd.concat([H,H_const])#H.con(H_const, ignore_index=True) #  %const. weak layer THICKNESS important for failure initiation criterion
    D_const =  pd.Series( [df_prof[rho].iloc[ind], np.sum(df_prof[rho].iloc[ind+1:-2] * df_prof[h].iloc[ind+1:-2] *  np.cos(alp+np.pi/180) ) \
                                             / np.sum( df_prof[h].iloc[ind+1:-2] *  np.cos(alp+np.pi/180))  ] )
    D = pd.concat([D,D_const])
    
    ### Slab layer prop
    E = Ela(D)
    TS = TSjj(D)

    ### 'average' slab properties
    rhoslab = np.sum(df_prof[h].iloc[0:ind]* df_prof[rho].iloc[0:ind]) \
                / np.sum(df_prof[h].iloc[0:ind])
    hslab = np.sum(df_prof[h].iloc[0:ind+1]) #%THICKNESS
    tau_g = rhoslab * 9.8 * hslab * np.sin( alp * np.pi/180)
    sigma_g = rhoslab*9.8 * hslab * np.cos( alp * np.pi/180)
    
    if calcFEM:
        print('get EBulk FEM-bulk not implemented') # get Ebulk line 89
    else:
        Eslab = np.sum(H[0:-2]*E[0:-2]) / np.sum(H[0:-2])
        # Monti's expression: skierloadE
    
    """WL properties"""
    rhowl = df_prof[rho].iloc[ind]
    
    if overwrite_scmod_strength:
        #print ('implement fu_tau_p_COJ15')
        tau_p , c_0 = fu_tau_p_CoJ15.fu_tau_p_CoJ15(index, df_P, WLopt = WLopt, pre_existing=pre_existing, dap_nr = dap_nr)
        #%% IMPLEMENT strengthening from fu_precstabindx.m ?
        print ('tau_p',tau_p,'c_0', c_0)
    else:
        
        tau_p = df_prof['snow shear strength (kPa)'].iloc[ind]*1000 # kPa in Pa
    
    wf = CERR(ini_rc, tau_p, Eslab) # %!rhoslab decreases when newsnow -> this increases wf (nonsense)
   
    #% wftemporal=CERR(ini_rc,tau_p,Ela(rhowl)); % if no FEM for bulk modulus, could use WL modulus MATLAB relikt
    
    """DAMAGE (natural release)"""
    ### calculate spontaneous release, if nap or pap/dap with new snow or rain
    ### Conway and Wilbur 1999, Snow Ava Form 2003, Capelli?, Reuter TCD
    ### burialdate='20141217'; dy='20150104'; %20160108
    if WLopt!='dap':
        if burialdate != pd.Timestamp('1970-01-01'):
            burial_stamp = pd.to_datetime(burialdate).replace(hour=0) # date of burial
            current_stamp = df_prof['date'].iloc[2] # todays stamp
            scmodstep = (df_met.timestamp[3]-  df_met.timestamp[2]).total_seconds()/3600 #interval of dataset in h
            
            tm = (current_stamp - burial_stamp).total_seconds()/3600 # timewindow in h

            tm = list(range(0,int(tm)+1))* int(scmodstep)
            
            
            mask_1 = (df_met.timestamp >= burial_stamp) 
            mask_2 =  (df_met.timestamp <= current_stamp)
            precrat = df_met.MS_Snow.loc[mask_1 & mask_2] # % snp: [smet. solid prec]=kg/m2/h
                
            ### Get natural instability due to precipitation
            ### 'buh' is spaceholder for not needed returned variables
            precstabMIN24,extm2failMIN24,tmcrit,buh,buh,buh,buh,buh,buh,buh = fu_precstabindx( burialdate, tm, precrat, scmodstep, alp)
        else:
            #% Stability metrics
            # precstabMIN24 = sig / precstr    # precipitation stability index  ConwayWilbur99 eq 6
            # extm2failMIN24 = (precstab-1) / np.append(np.NaN, np.abs(np.diff(precstabMIN24))) 
            # extm2failMIN24[extm2failMIN24 == 0] = np.NaN
            if WLopt == 'nap':
                height_opt = 'napup'
            elif WLopt == 'pap':
                height_opt = 'papup'
            precstabMIN24 = float(df_prof.loc[(df_prof['height_m'] == df_P.loc[index,height_opt]), 'Sn38' ].values[0]  )  # precipitation stability index  ConwayWilbur99 eq 6
            
            extm2failMIN24 = (precstabMIN24-1) / 0.05 # 0.05 = dSn/dt (per h) assuming small change of 0.05 within one hour 
            tmcrit = np.nan
    else:
        precstabMIN24  = np.nan
        extm2failMIN24 = np.nan
        tmcrit         = np.nan

    ### Natural stability index Sn
    Sn = tau_p / tau_g

    """INITIATION"""
    ### Calc initiation criterion, provide data for slab, weak layer and basal layer
    msswl,skierloadE,equiHslab,mssANA = pro_rb15.get_S_rb15_v2(H,D,E,tau_p,alp,1,calcFEM)
    Sk_ana = tau_p / (tau_g + mssANA);  #% McClung & Sz 99, Monti etal 16
    Sk_fem = tau_p / (tau_g + msswl) ;   #% Rb & Sz 18 GRL append. 

    ### Skier stability index (Sk38) at a slope (estimate deltatauxy, otherwise use get_S_rb15.m)
    angle = [30,35,38]
    angle_ind = []
    for element in angle:
        x = element == alp
        angle_ind.append(x)
    phi=np.array([55.49,54.77,54.34])
    angle_ind = np.array(angle_ind)
    angle = np.max(phi*angle_ind)
    
    if np.sum(angle_ind)==1:
            deltau = 2 * 500 * np.cos(angle*np.pi/180) *np.sin(angle*np.pi/180)**2 \
                     * np.sin((angle+alp)*np.pi/180) / np.pi  / hslab
    else:
        deltau=150 / (hslab / np.cos(alp*np.pi()/180)) # short for alp 35 and lineload 500N/m (Schweizer internal report says H is thickness, but it is height, so correct with cos!)

    Sk = tau_p / (tau_g + deltau)#; % Sz internal repor

    """PROPAGATION"""
    ### Calc critical crack length
    ac_vh16,wf,ac_si06,ac_ga17 = pro_vh16.get_ac_vh16_v2(Eslab, rhoslab, hslab, wf, alp, rho_wl=rhowl,tau_p=tau_p)


    """CHECK FOR SLAB SUPPORT AND SUSTAINED (DYNAMIC) PROPAGATION"""
    ### Calc tensile support of the slab
    if calcFEM == 1:
        print('[I]  FEM batch not implemented')
        #  [tcr,tenscrit,slablaythick]=get_ts_rb17(H(1:end-2),D,E,TS,Ltot,alp,rc); % check threshold for critical in the routine
    
    ### Analytic solution 
    ### Is the fracture speed an indicator for sustained propagation?

    """Prepare output of fu_instab()"""
    DAM=np.array([Sn,precstabMIN24,extm2failMIN24,tmcrit])
    INI = np.array([tau_p, c_0, Sk, Sk_ana, mssANA, Sk_fem, msswl,sigma_g])
    PRO=np.array([ac_si06, ac_vh16, wf, ac_ga17])
    #PRO=[ac_vh16, wf];
    DYN=tcr
    # disp(['INI S: ',num2str(tau_p/mssANA)])
    # disp(['PRO ac: ',num2str(ac_vh16)])
    return DAM,INI, DYN, PRO