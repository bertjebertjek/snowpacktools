import numpy as np
import pandas as pd

def fu_tau_p_CoJ15(index, df_P,opt = 'per',WLopt = 'pap', pre_existing = False, dap_nr = False):
    ''' Function to to calcuate the evolution of weak layer shear strength according to the 
    strengthening observed in faceted WLs (Conlan & Jamieson 2015, CGJ)

    Input:
        in df:
            dy: timestamp
            dynoload: days without addtional added load on the snowpack (days without snowfall)
            ini_tau_p: yesterday's shear strength or initial shear strength (if WL is less than a day old)
            c_0: base rate for shear strength growth law (only needed if WLage > 1 )
            WLage in days
            rhowl: weak layer density
        opt: persistent or non-persistent weak layers
            case1 : 'per'
            case2 : 'nonper'
    Output:
        tau_p: initial shear strength
        c_0: base rate for shear strength growth law (only needed if WLage > 1 )
    '''

    #%% define strengthening exponent
    #% other resources: DF: Reuter TCD (100/dy), [vH Miller13 (.11) !densific.!]; RG: vH Miller13 (.2),
    #% Milleretal03 (.18); SH: [Hortonetal14 (.26) !infrequent!]; DH: [vHMiller13 (.21) !sifted!]
    
    # handling different WL options
    if opt == 'per':
        strexp = 0.26
        ini_tau_p = 200  # ini_tau_p=SSjjper(rhowl);
    elif opt == 'nonper':
        strexp = 0.26 
        ini_tau_p = 200 #%ini_tau_p=SSjjnon(rhowl);
    else:
        print('per/nonper option needed as input')

    if WLopt == 'pap':
        WLup = 'papup'
        burialdate = df_P['papburial'][index] # 'papburial'
    elif WLopt == 'nap':
        WLup = 'napup'
        burialdate =  df_P['napburial'][index] #'napburial'
        
    elif WLopt == 'dap':
        WLup = df_P['dapup'][index][dap_nr]
        burialdate = df_P['dapburial'][index][dap_nr] #
        ini_tau_p = df_P['dapINI_tau_p'][index][dap_nr]
        if ini_tau_p == None:
                ini_tau_p = 0
        c_0 = df_P['dapINI_c_0'][index][dap_nr]

    if pre_existing:
        if WLopt == 'pap':
            ini_tau_p = df_P['papINI_tau_p'][index-1]
            if ini_tau_p == None:
                ini_tau_p = 0
            c_0 = df_P['papINI_c_0'][index-1]

        elif WLopt == 'nap':
            ini_tau_p = df_P['napINI_tau_p'][index-1]
            if ini_tau_p == None:
                ini_tau_p = 0
            c_0 = df_P['napINI_c_0'][index-1]

    
    # laws
    def SScj(dys, c_0):
        SScj = c_0*(dys) ** strexp #exp. strengthening law
        return SScj
        
    #SScj=@(dys,c_0)(c_0.*(dys).^strexp);    %exp. strengthening law
    #SSjjper=@(rho)(18.5e3*(rho/917).^2.11); % shear strength persistent grains (JamJohn2001)
    #SSjjnon=@(rho)(14.5e3*(rho/917).^1.73); % shear strength non-persistent grains (JamJohn2001)

    if burialdate != pd.Timestamp('1970-01-01'):
        WLage = (df_P.dy[index] - burialdate.date()).days   # duration of burial in days
    else:
        WLage = 2
    #%% initial strength, if initial values (ini_tau_p, c_0) not provided
    #% DF/SH/siftedFC: 1kPa (Reuteretal19 TCD)
    if WLage < 1:
        tau_p = ini_tau_p         # initial shear strength
        c_0 = tau_p*2**(-strexp)  # base rate of the strengthening law
    else:
        if df_P['dynoload'][index] == 0:
            #% new snow loading, i.e. strength increases
            arg = 2
            c_0 = ini_tau_p
        else:
            #% no new snow, strength increases more gradually
            #% don't change c_0, the exponential law will do the job
            arg = 2 + df_P.dynoload[index]
        
        if 'c_0' not in locals():
            c_0 = ini_tau_p

        tau_p = SScj(arg, c_0) #[Pa] use Conlan Jamieson 2015
    
    if ini_tau_p == None:
        ini_tau_p = 0   
        
          
    return tau_p, c_0


              