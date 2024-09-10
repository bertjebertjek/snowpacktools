################################################################################
# Copyright 2022 Avalanche Warning Service Tyrol                               #
################################################################################
# This is free software you can redistribute/modify under the terms of the     #
# GNU Lesser General Public License 3 or later: http://www.gnu.org/licenses    #
################################################################################

import warnings
import pandas as pd
import numpy as np

def assign_aps(df_P, config):
    """Assign avalanche problems (APs) based on collected WLs"""
    
    """Thresholds"""
    debug       = int(config.get('AVAPRO', 'debug'))
    dysisomax   = int(config.get('AVAPRO-THOLDS-WL', 'dysisomax'))

    if config.get('AVAPRO', 'scmopt') == 'snp':
        thold_scmopt = 'AVAPRO-THOLDS-APS-SNP'
    elif config.get('AVAPRO', 'scmopt') == 'cro':
        thold_scmopt = 'AVAPRO-THOLDS-APS-CRO'
    else:
        raise ValueError('[E]   Check scmopt in ini file, only snp or cro valid. your input is:', config.get('AVAPRO', 'scmopt'))

    inithrshnap     = int(config.get(thold_scmopt, 'inithrshnap'))     # do not constrain
    propthrshnap    = int(config.get(thold_scmopt, 'propthrshnap'))
    damthrshnap     = int(config.get(thold_scmopt, 'damthrshnap'))     # snp: .32m
    damthrshpap     = int(config.get(thold_scmopt, 'damthrshpap'))     # snp:18 8 18 hours
    # inithrshdap   = 2x inithrshpap; propthrshdap= 2x propthrshpap;
    inithrshpap     = float(config.get(thold_scmopt, 'inithrshpap'))   # snp:1;
    propthrshpap    = float(config.get(thold_scmopt, 'propthrshpap'))  # snp:.24 .318 .4;
    drftthrsh       = int(config.get(thold_scmopt, 'drftthrsh'))       # snp.5[m]; min height of snow drifts


    """Look at PAPs of season"""
    df_P['papex_sele'] = df_P['papex']*df_P['papcalc']
    df_P.loc[df_P['papex_sele']!=1, ['papex_sele'] ] = 0
    df_P['pap_isrel'] = [np.nan for _ in range(len(df_P))]
    df_P['papDAM_precstabMIN24_isrel'] = [np.nan for _ in range(len(df_P))]
    df_P['papDAM_extm2failMIN24_isrel'] = [np.nan  for _ in range(len(df_P))]
    df_P['papINI_isrel'] = [np.nan for _ in range(len(df_P))]
    df_P['papPRO_isrel'] = [np.nan for _ in range(len(df_P))]

    for index,row in df_P.iterrows():
        if df_P.loc[index,'papex_sele'] == 1:
            if (df_P.loc[index,'papDAM_precstabMIN24'] > 0) & (df_P.loc[index,'papDAM_extm2failMIN24'] >0):
                df_P.loc[index,'papDAM_precstabMIN24_isrel'] = df_P.loc[index,'papDAM_precstabMIN24'] 
                df_P.loc[index,'papDAM_extm2failMIN24_isrel']  = df_P.loc[index,'papDAM_extm2failMIN24'] 
                df_P.loc[index,'papINI_isrel']  = df_P.loc[index,'papINI_tau_p']  / df_P.loc[index,'papINI_mssANA'] 
                df_P.loc[index,'papPRO_isrel']  = df_P.loc[index,'papPRO_ac_vh16'] 
    ### ??papPRO_is_rel len 3 kürzer (1. nicht gecatched; matlab 2 nans ok)
    df_P['papdycrit'] = [0 for _ in range(len(df_P))]     
    
    ### Natural release criteria
    df_P['papex_sele_natural'] = np.where ((df_P['papex_sele'] == 1) & ((df_P['papDAM_extm2failMIN24_isrel'] > damthrshpap) | (df_P['papPRO_isrel'] >propthrshpap)), 0, df_P['papex_sele'] )
    df_P['papdycrit_natural'] =  np.where ((df_P['papex_sele'] == 1) & ((df_P['papDAM_extm2failMIN24_isrel'] <= damthrshpap) & (df_P['papPRO_isrel'][0]<= propthrshpap)), 1,0 )    
    
    ### Artificial (triggered) release criteria
    df_P['papex_sele_trigger'] = np.where ((df_P['papex_sele'] == 1) & ((df_P['papPRO_isrel'] > propthrshpap) | (df_P['papINI_isrel'] >  inithrshpap)  ), 0, df_P['papex_sele']  )
    df_P['papdycrit_trigger'] =  np.where ((df_P['papex_sele'] == 1) & ((df_P['papPRO_isrel'] <= propthrshpap)  & (df_P['papINI_isrel'] <= inithrshpap)  ),1,0)    

    print('[i]  PAPs (trigger):',np.sum(df_P['papex_sele_trigger']))
    print('[i]  PAPs (natural):',np.sum(df_P['papex_sele_natural']))


    """Look at NAPs of season"""
    df_P['napex_sele'] = [np.nan for _ in range(len(df_P))]
    df_P['napDAM_extm2failMIN24_isrel'] = [np.nan  for _ in range(len(df_P))]
    df_P['napINI_isrel'] = [np.nan for _ in range(len(df_P))]
    df_P['napPRO_isrel'] = [np.nan for _ in range(len(df_P))]
    df_P['napDAM_precstabMIN24_isrel'] = [np.nan for _ in range(len(df_P))]

    for index,row in df_P.iterrows():  
            df_P.loc[index,'napex_sele']  = df_P.loc[index,'napex'] * df_P.loc[index,'napcalc'] 

            if (~np.isnan(df_P['napex_sele'][index])) & (df_P['napex_sele'][index] != 0) : 
                    df_P.loc[index,'napex_sele'] = 1
                    df_P.loc[index,'napDAM_precstabMIN24_isrel'] = df_P.loc[index,'napDAM_precstabMIN24']
                    df_P.loc[index,'napDAM_extm2failMIN24_isrel'] = df_P.loc[index,'napDAM_extm2failMIN24']
                    df_P.loc[index,'napINI_isrel']  = (df_P.loc[index,'napINI_tau_p'] / df_P.loc[index,'napINI_mssANA'])
                    df_P.loc[index,'napPRO_isrel'] = df_P.loc[index,'napPRO_ac_vh16'] 
                    if np.iscomplexobj(df_P.loc[index,'napPRO_isrel']): 
                        df_P.loc[index,'napPRO_isrel']  = np.nan

    ### Natural release criteria
    df_P['napex_sele_natural'] = np.where ((df_P['napex_sele'] == 1) & ((df_P['napDAM_extm2failMIN24_isrel'] > damthrshnap) | (df_P['napPRO_isrel'] >propthrshnap  )), 0, df_P['napex_sele'] )  
    df_P['napdycrit_natural'] =  np.where ( (df_P['napex_sele']== 1) & ((df_P['napDAM_extm2failMIN24_isrel'] <=  damthrshnap) & (df_P['napPRO_isrel'] <= propthrshnap)), 1,0 )  

    ### Artifical (triggered) release criteria
    df_P['napex_sele_trigger'] = np.where ( (df_P['napex_sele'] == 1) & ((df_P['napINI_isrel'] > inithrshnap) | (df_P['napPRO_isrel'] >  propthrshnap)), 0, df_P['napex_sele']  )  
    df_P['napdycrit_trigger'] =  np.where ( (df_P['napex_sele'] == 1) & ((df_P['napINI_isrel'] <= inithrshnap) & (df_P['napPRO_isrel'] <= propthrshnap)),1,0)  
    # print(np.sum(df_P['napex_sele'] == 1))
    # print(np.sum((df_P['napINI_isrel'] <= inithrshnap) & (df_P['napPRO_isrel'] <= propthrshnap)))
    # print(np.sum(df_P['napdycrit']))

    print('[i]  NAPs (trigger):',np.sum(df_P['napex_sele_trigger']))
    print('[i]  NAPs (natural):',np.sum(df_P['napex_sele_natural']))

    """Look at DAPs of season"""
    df_P['dap_isrel'] = [[0]*50 for _ in range(len(df_P))]
    df_P['dapDAM_precstabMIN24_isrel'] = [[np.nan]*50 for _ in range(len(df_P))]
    df_P['dapDAM_extm2failMIN24_isrel'] = [[np.nan]*50  for _ in range(len(df_P))]
    df_P['dapINI_isrel'] = [[np.nan]*50 for _ in range(len(df_P))]
    df_P['dapPRO_isrel'] = [[np.nan]*50 for _ in range(len(df_P))]
    df_P['dapup_isrel'] =  [[np.nan]*50 for _ in range(len(df_P))]
    df_P['dapex_isrel'] = [[np.nan]*50 for _ in range(len(df_P))]
    df_P['dapex_isrel_red'] = [0 for _ in range(len(df_P))]
    
    for index,row in df_P.iterrows():
        # if debug: print('[D]    DAP exists: ', df_P['dapex'][index])

        if ~np.isnan(df_P['dapex'][index]):
            ### No. of DAPs of current/index day
            nr_daps = len(df_P['dapINI_tau_p'][index])
            
            for ele in range(nr_daps):
                if ~np.isnan(df_P['dapup'][index][ele]):
                    df_P.loc[index, 'dapex_isrel'][ele] = 1
                    df_P.loc[index,'dapINI_isrel'][ele] = df_P.loc[index,'dapINI_tau_p'][ele] / df_P.loc[index,'dapINI_mssANA'][ele] # %ini_ana
                    df_P.loc[index,'dapPRO_isrel'][ele] = df_P.loc[index,'dapPRO_ac_vh16'][ele] # ac_vh16
                    
                    ### Set 0 to np.nan
                    if df_P.loc[index,'dapINI_isrel'][ele] == 0:
                        df_P.loc[index,'dapINI_isrel'][ele] = np.nan
                    if df_P.loc[index,'dapPRO_isrel'][ele] == 0:
                        df_P.loc[index,'dapPRO_isrel'][ele] = np.nan
        
            # print(np.sum( ~np.isnan(df_P['dapex_isrel'][index])))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                df_P.loc[index,'dapINI_isrel'] = np.nanmin(df_P.loc[index,'dapINI_isrel'])
                df_P.loc[index,'dapPRO_isrel'] = np.nanmin(df_P.loc[index,'dapPRO_isrel'])
            df_P.loc[index,'dapex_isrel_red'] = int(np.sum( ~np.isnan(df_P.loc[index,'dapex_isrel']))>0)


    """DAP is secondary release of new snow, wind slab AND/OR persistent (PAP)"""
    ### Natural release criteria
    df_P['dapex_sele_natural'] = np.where((df_P['dapex_isrel_red']== 1) & (~ (df_P['napdycrit_natural'] == 1)),  0, df_P['dapex_isrel_red'])
    
    ### Artifical (triggered) release criteria
    df_P['dapex_sele_trigger'] = np.where((df_P['dapex_isrel_red']== 1) & (~((df_P['papdycrit_trigger']  == 1) | (df_P['napdycrit_trigger'] == 1))),  0, df_P['dapex_isrel_red'])

    print('[i]  DAPs (trigger):', np.sum(df_P['dapex_sele_trigger']))
    print('[i]  DAPs (natural):', np.sum(df_P['dapex_sele_natural']))
    
    """Look at WAPs of season"""
    ### Define a variable that indicates if LWC increased with respect to the previous day
    # df_P.loc[df_P['wapLWC'].isna(), 'wapLWC'] = 0 # (format for plot)
    i = np.arange(0,len(df_P['wapLWC'])-1)
    iplus = np.arange(1,len(df_P['wapLWC']))
    df_P['incr'] = df_P['wapLWC'].copy()
    df_P.loc[iplus, 'incr'] = np.where(((df_P['wapLWC'].values[iplus]-df_P['wapLWC'].values[i])>0.0001), 1, np.nan)
    df_P['incr'].values[0] = np.nan
    
    ### Count days in a row the snowpack is isothermal based on defining isothermal state in find_aps()
    df_P.loc[df_P['wapISO'].isna(), 'wapISO'] = 0 
    df_P['dysio'] = [np.nan for _ in range(len(df_P))]

    df_P.loc[0,'dysio'] = 0
    for ii in range(1,len(df_P)):
        if df_P.loc[ii,'wapISO']:
            df_P.loc[ii,'dysio'] = df_P.loc[ii-1,'dysio'] + 1
        else:
            df_P.loc[ii,'dysio'] = 0

    ### Format 'waponset'
    df_P.loc[df_P['waponset'].isna(), 'waponset'] = pd.NaT
    
    # df_P['wapex_sele'] = np.where( (df_P['dy'] >= df_P['waponset']) &  df_P['incr'] & (df_P['dysio'] <dysisomax ), 1, np.nan)
    df_P['wapex_sele'] = np.where( (pd.to_datetime(df_P['dy']) >= df_P['waponset'])  & (df_P['dysio'] <= dysisomax), 1, np.nan)
    print('[i]  WAPs:', np.sum(df_P['wapex_sele']))

    # df_P['wapLWC_isrel'] = np.where(~df_P['wapex_sele'].isna(), df_P['wapLWC'], np.nan )
    # df_P['wapLWC_isrel'] = np.where(df_P['wapcycle'] ==1, df_P['wapLWC_isrel'] / lwcthrsh_0, df_P['wapLWC_isrel'])
    # df_P['wapLWC_isrel'] = np.where(df_P['wapcycle'] > 1, df_P['wapLWC_isrel'] / lwcthrsh_1, df_P['wapLWC_isrel'])
    # print('[i]  WAPs:',len(df_P['wapLWC_isrel'][df_P['wapLWC_isrel'] > 0]))
    

    """Wind (WSAP/winex) based on VW and loose snow"""
    # df_P['winex_sele_trigger'] = df_P['winex']
    print('[i]  WSAPs:', np.sum(df_P['winex']))

    ### (WE ALSO NEED TO GET THE WIND INDICATOR FOR WIND LOADING (rather than new snow loading -- pure wind slab case))
    ### (Days with neither pap nor nap, but wind loading)
    #? sele=find(P.napwindonly)
    #? df_P['winex'] = df_P['napwindonly']
    #? df_P['winex'] = np.nan
    #? df_P['winex'] = np.where ( (df_P['drft']> drftthrsh ) , 1, 0 )

    return df_P