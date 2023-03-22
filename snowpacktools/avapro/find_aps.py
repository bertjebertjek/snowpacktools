
################################################################################
# Copyright 2022 Avalanche Warning Service Tyrol                               #
################################################################################
# This is free software you can redistribute/modify under the terms of the     #
# GNU Lesser General Public License 3 or later: http://www.gnu.org/licenses    #
################################################################################

# import warnings
# warnings.filterwarnings('ignore')

import pickle
import numpy as np
import pandas as pd
import datetime

from snowpacktools.avapro import fu_instab as pro_instab
from snowpacktools.avapro import fu_tau_p_CoJ15 as fu_tau_p_CoJ15
from snowpacktools.snowpro import snowpro

def find_aps(config, pro_path, smet_path):
    """Search for weak layers and identify the corresponding avalanche problem (AP) or AP type
    Input:
        pro_path:
        smet_path:
        
    Output:
        nap:            Non-persistent AP, keep for +1day unless HN24 exceeds threshold
        pap:            Persistent AP, keep, becomes dap if more recent pap appears
        dap (yymmdd):   Deep (persistent) AP, was pap before, has label of WL burial date
        
    """

    """Parameters for research applications"""
    read_in = 1
    dump_profs_and_smet = 0 # Needs to be set to 1 once before setting read_in = 0

    ### initialization_type: For simulations initialized with profiles, a startup routine based on RTA is used Can be station or profile.
    initilization_type    = config.get('AVAPRO', 'initilization_type')
    season_start = config.get('AVAPRO', 'season_start') 
    season_end   = config.get('AVAPRO', 'season_end')
    drytime = int(config.get('AVAPRO', 'drytime'))
    wettime = int(config.get('AVAPRO', 'wettime'))

    """Thresholds"""
    drftthrsh   = float(config.get('AVAPRO-THOLDS-WL', 'drftthrsh'))     # min height of snow drifts [m], snp=:0.4m, cro:.2
    nsthrsh     = float(config.get('AVAPRO-THOLDS-WL', 'nsthrsh'))       # def.: .2,amount of new snow [m] (within 48 or 24 hr), to keep a non-persis. problem for one more day   
    # calcFEM     = int(config.get('AVAPRO-THOLDS-WL', 'calcFEM'))         # =1; means use FEM model to calcuate shear stress at WL, slab modulus, slab tensile stresses
    # owSCtaup    = int(config.get('AVAPRO-THOLDS-WL', 'owSCtaup'))        # =1; means overwrite snow cover shear strength (tau_p)
    # lookat      = int(config.get('AVAPRO-THOLDS-WL', 'lookat'))          # Say if you want to see a graph    
    # minfracswetsnow=0.7                                                # min fraction of snow height for wap, def.: 0.7 * max snow height during season
    # thresholds_aps['release'] = release

    if read_in == 1:
        print('[I]  Generate list_df_pro (list of profiles)')
        list_df_pro, meta_dict = snowpro.read_pro(pro_path)
        print('[I]  Generate df_met')
        df_met = snowpro.get_smet_df(smet_path)
        print('[I]  Slope Angle: {}, Aspect: {}'.format(meta_dict['SlopeAngle'],meta_dict['SlopeAzi']))

        # ---- DUMPING PKL FILES OF PROFILES ---- #
        if dump_profs_and_smet:
            with open('./list_df_pro', 'wb') as f:
                pickle.dump(list_df_pro, f)

            with open('./df_met', 'wb') as f:
                pickle.dump(df_met, f)
    else :
        # ---- READING PKL FILES OF PROFILES ---- #
        print('[I]  Read in existing list_df_pro and df_met')
        list_df_pro = pickle.load(open('./list_df_pro', "rb"))
        #meta_dict = pickle.load(open('./meta_dict', "rb"))
        df_met = pickle.load(open('./df_met',"rb"))

    # ---- Reduce met file to season ---- #
    df_met['timestamp'] = pd.to_datetime(df_met['timestamp'])
    
    if initilization_type=='station': 
        mask_season = (df_met['timestamp'] > season_start) & (df_met['timestamp'] < season_end)
        df_met_red= df_met.loc[mask_season].reset_index(drop=True)
    else:
        df_met_red=df_met
    
    # ---- Reduce profiles and met file to dry and wet time ---- #
    mask_dry = (df_met_red['timestamp'].dt.hour == drytime ) 
    mask_wet = (df_met_red['timestamp'].dt.hour == wettime )
    df_met_red = df_met_red.loc[mask_dry | mask_wet]
    df_met_red = df_met_red.reset_index(drop=True)

    season_list = []
    season_list_red = []

    if initilization_type=='station':
        ### Reduce profile list to season start/end
        for df in list_df_pro:
            mask = (df.date > season_start) & (df.date  < season_end)
            if mask[0]:
                season_list.append(df)
    else:
        for df in list_df_pro:
            season_list.append(df)

    ### Reduce season_list to dry&wet time:
    for df in season_list:
        ### Correct profile with slope angle for a consistent calculation of Sn38, Sk38, ...
        df['height_m']    = df['height_m'] * np.cos(float(meta_dict['SlopeAngle'])*np.pi/180)
        df['thickness_m'] = df['thickness_m'] * np.cos(float(meta_dict['SlopeAngle'])*np.pi/180)
        df['bottom']      = df['bottom'] * np.cos(float(meta_dict['SlopeAngle'])*np.pi/180)

        if df.date.iloc[0].hour == drytime or df.date.iloc[0].hour == wettime:
            season_list_red.append(df)
            
    P_varibales = [
        'dy',  #time to look up the morning profile,
        'dynoload', # day since last significant (define threshold) snowfall 
        'napex', 'papex','dapex', #indicator for (non-)persistent
        'napup', 'papup','dapup', #upper end of layer
        'naplo', 'paplo', 'daplo', #lower end of layer
        'napSLdep', 'papSLdep','dapSLdep', #%slab depth
        'napSLrho', 'papSLrho','dapSLrho', #slab density
        'napWLrho', 'papWLrho','dapWLrho', #%WL density
        'napcalc','papcalc',    #decision to calc instability
        'napgt', 'papgt','dapgt', #WL grain type
        'napburial', 'papburial','dapburial', # WL burial date
        'napwindonly', #pure wind slab problem, napex is set 2!
        ### WAP
        'wapex', 'wapcalc', 'wapLWC','wapSWE','wapTWAT', 
        'wapISO', 'wapWLdry','waponset','wapcycle',
        ### NAP 
        'napDAM_Sn','napDAM_precstabMIN24','napDAM_extm2failMIN24',
        'napDAM_tmcrit','napINI_tau_p','napINI_c_0','napINI_Sk','napINI_Sk_ana',
        'napINI_mssANA','napINI_Sk_fem', 'napINI_msswl','napINI_sigma_g',
        'napDYN','napPRO_ac_si06','napPRO_ac_vh16','napPRO_wf', 'napPRO_ac_ga17',
        ### PAP
        'papDAM_Sn','papDAM_precstabMIN24','papDAM_extm2failMIN24',
        'papDAM_tmcrit','papINI_tau_p','papINI_c_0','papINI_Sk','papINI_Sk_ana',
        'papINI_mssANA','papINI_Sk_fem', 'papINI_msswl','papINI_sigma_g',
        'papDYN','papPRO_ac_si06','papPRO_ac_vh16','papPRO_wf', 'papPRO_ac_ga17',
        ### DAP
        'dapDAM_Sn','dapDAM_precstabMIN24','dapDAM_extm2failMIN24',
        'dapDAM_tmcrit','dapINI_tau_p','dapINI_c_0','dapINI_Sk','dapINI_Sk_ana',
        'dapINI_mssANA','dapINI_Sk_fem', 'dapINI_msswl','dapINI_sigma_g',
        'dapDYN','dapPRO_ac_si06','dapPRO_ac_vh16','dapPRO_wf', 'dapPRO_ac_ga17',
        ### WSAPs
        'winex','winex_count',
        ### Meteo-Data
        'hn24','hn48','drft','hs',
        'snowclim']


    """Create empty dataframe df_P"""
    if df_met_red.loc[0,'timestamp'].hour == 15:
        df_met_red = df_met_red.iloc[1:,:]
        season_list_red = season_list_red[1:]
        df_met_red = df_met_red.reset_index(drop=True)

    date = df_met_red['timestamp'][::2]
    date = pd.to_datetime(date).apply(lambda x: x.date())
    
    ### Generate empty df with len of date
    if initilization_type=='station':
        df_P  = pd.DataFrame(np.nan, index= range(len(date)), columns=P_varibales)
        df_P ['dy'] = date.values   ### add timestamp from metfile to final df
    else:
        ### Add one entry at the beginning for profile initiliazation!
        df_P  = pd.DataFrame(np.nan, index= range(len(date)+1), columns=P_varibales)
        df_P.loc[0:,'dy'] = date.values[0]-pd.Timedelta(days=1)
        df_P.loc[1:,'dy'] = date.values

    df_P['dapup'] = [[] for _ in range(len(df_P))]
    df_P['daplo'] = [[] for _ in range(len(df_P))]
    df_P['dapSLdep'] = [[] for _ in range(len(df_P))]
    df_P['dapSLrho'] = [[] for _ in range(len(df_P))]
    df_P['dapWLrho'] = [[] for _ in range(len(df_P))]
    df_P['dapburial'] = [[] for _ in range(len(df_P))]
    df_P['dapgt'] = [[] for _ in range(len(df_P))]
    df_P['dapDAM_Sn'] = [[] for _ in range(len(df_P))]
    df_P['dapDAM_precstabMIN24'] = [[] for _ in range(len(df_P))]
    df_P['dapDAM_extm2failMIN24'] = [[] for _ in range(len(df_P))]
    df_P['dapDAM_tmcrit'] = [[] for _ in range(len(df_P))]
    df_P['dapINI_tau_p'] = [[] for _ in range(len(df_P))]
    df_P['dapINI_c_0'] = [[] for _ in range(len(df_P))]
    df_P['dapINI_Sk'] = [[] for _ in range(len(df_P))]
    df_P['dapINI_Sk_ana'] = [[] for _ in range(len(df_P))]
    df_P['dapINI_mssANA'] = [[] for _ in range(len(df_P))]
    df_P['dapINI_Sk_fem'] = [[] for _ in range(len(df_P))]                                
    df_P['dapINI_msswl'] = [[] for _ in range(len(df_P))]
    df_P['dapINI_sigma_g'] = [[] for _ in range(len(df_P))]
    df_P['dapPRO_ac_si06'] = [[] for _ in range(len(df_P))]
    df_P['dapPRO_ac_vh16'] = [[] for _ in range(len(df_P))]
    df_P['dapPRO_wf'] = [[] for _ in range(len(df_P))]
    df_P['dapPRO_ac_ga17'] = [[] for _ in range(len(df_P))]
    df_P['dapDYN'] = [[] for _ in range(len(df_P))]
        
    df_P['napgt'] =  [[np.nan, np.nan, np.nan] for _ in range(len(df_P))]
    df_P['papgt'] =  [[np.nan, np.nan, np.nan] for _ in range(len(df_P))]

    cycle_nr = 1 ### Tracking of wet cycles
    runno = 1

    """Initialization"""
    if initilization_type=='profile':
        print("[I]  Startup-Algorithm: Initializing PAPs and DAPs with RTA")
        df_P = initialize_with_RTA(config,season_list,df_met,df_P)

    """Iteration"""
    print("[I]  Start of iteration")
    # for index, row in islice(df_P.iterrows(), 1, None):
    for index in range(1,len(df_P)):
        ### Index handling
        if initilization_type=='station':
            ind_dry = index*2 # drytime = 6
            ind_wet = ind_dry +1
        else:
            ind_dry = (index-1)*2 # drytime = 6
            ind_wet = ind_dry +1

        ### dy is the date
        dy = df_met_red.timestamp.iloc[ind_dry]
        ### Address snowheight from smet to final df acording to date
        df_P.loc[index, 'hn24'] = df_met.loc[(df_met["timestamp"] == dy), 'HN24'].values[0] 
        ### [m] height of snow drifts within the past 24h before 8:00 a.m. ???6 am ?
        df_P.loc[index, 'drft'] = df_met_red.loc[ind_dry, 'wind_trans24']
        ### Sum up 24h snowfall to 48h
        df_P.loc[index , 'hn48']  = np.sum( [df_P.hn24.iloc[index] ,df_P.hn24.iloc[index-1]])

        ### Define start value
        if runno == 1:  
            hsprotdy = df_met_red.HS_mod[ind_dry].copy()
        
        ### Days without adding load / days since last snowfall (preparation for strength calculation AND for natural release) 
        if  df_P.loc[index, 'hn24'] > nsthrsh or df_P.loc[index, 'drft'] > drftthrsh: 
            df_P.loc[index, 'dynoload'] = 0 #  signif. loading -> strength remains (see tau_p module)
        else:
            df_P.loc[index, 'dynoload'] = df_P.loc[index-1, 'dynoload'] + 1 #  % no signif. snowfall -> strength increases

        ### Read-in of drytime profile  
        df_prof = season_list_red[ind_dry]
        
        if df_prof.loc[0,'date']!= dy:
            print('date of season_list_red (profiles) -',df_prof.date[0],'- does not match date of df_P-',dy)
            break
        
        ### var 'po' refers to height of this layer (height [> 0: top, < 0: bottom of elem.] (m)'# just for better readability)
        po= 'height_m'

        if  df_prof[po][0] > 2 * nsthrsh:     
            ### minimum height of snowpack to start procedure
            ## hsprotdy = df_met_red.HS_mod[0]
            hsproydy = hsprotdy
            hsprotdy = df_prof[po][0]
            deltahs = hsprotdy - hsproydy
            ### Build shabby snow transport indicator
            
            """Find prior problems/WLs in today's profile"""
            ### matlab line 143
            df_P = find_prior_problems(config,index,df_prof,df_met,df_P)
            
            """Identify new problems/WLs in profile"""
            ### Mainly done in postprocessing - probably only one part needed in the end
            if ind_wet ==len(df_met_red.MS_Water): # catching boundary value problem
                pass
            else:
                df_P, cycle_nr = identify_wet_ap(config,index,ind_wet,season_list_red,df_met_red,df_P,cycle_nr)  

                """New NAP or PAP"""
                ### Check new snow and the prior surface
                ### We only enter if NO prior non-persistent problem pre-exists
                
                ### COMMENT ON WIND SLAB: This could also work for wind slabs covering WLs, if we use deltahs instead of hn24
                ### but snowpack flatfield sim do not in-/decrease snow height when 24h_wind_drift>0 !?          
                df_P = identify_new_nap_or_pap(config,index,ind_dry,season_list_red,df_met_red,df_met,df_P,deltahs)

        runno = runno+1
    return df_met_red, df_P, meta_dict


def find_lay_prof(config,df_prof,Lup,Lgt,wlopt):
    """Find vertical location of WL of prior profile in current profile
    - Make sure class (grain types) match ('per' vs 'nonper')
    """

    """Parameters"""
    debug       = int(config.get('AVAPRO', 'debug'))

    """Parameter initalization"""
    Lupnew = np.NaN
    Llonew = np.NaN
    Lgtnew = ['np.NaN']
    Lrho = np.NaN
    SLdep = np.NaN
    SLrho = np.NaN

    po = 'height_m'
    idx = sum((df_prof[po]-Lup)>0) +1 # index of WL
    selgt = df_prof.graintype[idx-1:idx+4] # graintypes in layers around WL
    selgt_comp = [] # check if prior graintype still in prof
    for ele in selgt:
        if ele[0] == Lgt[0]:
            selgt_comp.append(True)
        else:
            selgt_comp.append(False)

    
    if sum(selgt_comp)>0:
        ### grain type is still the same --> find straight away
        Lupnew = df_prof[po][idx+selgt_comp.index(True)-1]
        if Lupnew > df_prof[po][-1:].values:
            Llonew = df_prof[po][idx  +selgt_comp.index(True)]
            Lgtnew = Lgt
            SLdep = df_prof[po][0]-Lupnew        
            SLrho = np.sum(df_prof['density'][0:idx+selgt_comp.index(True)]  \
                            * df_prof['thickness_m'][0:idx+selgt_comp.index(True)]) \
                    / np.sum(df_prof['thickness_m'][0:idx+selgt_comp.index(True)])
            Lrho = df_prof['density'][idx+selgt_comp.index(True)-1]
            layfound = 1
        else:
            layfound = 0
    else:
        ### grain type changed, so look for other (non-)peristent types
        if wlopt == 'per':
            #print('check here per', index)
            id_p = selgt.str.contains('SH', regex=False) + selgt.str.contains('FC',regex = False) \
                        + selgt.str.contains('FCxr', regex=False) + selgt.str.contains('DH',regex =False)
            
            id_p = id_p.reset_index(drop=True)
            
            if sum(id_p.values)>0:
                if debug: print('[D]    WL grain type changed, looking for other persistent GT')
                Lupnew = df_prof[po][idx + np.where(id_p)[0][0]]
                if Lupnew > df_prof[po][-1:].values:
                    Llonew = df_prof[po][idx + np.where(id_p)[0][0]] 
                    Lgtnew = df_prof['graintype'][idx + np.where(id_p)[0][0]] 
                    SLdep = df_prof[po][0]-Lupnew        
                    SLrho = np.sum(df_prof['density'][0:idx+np.where(id_p)[0][0]]  \
                                    * df_prof['thickness_m'][0:idx+np.where(id_p)[0][0]]) \
                            / np.sum(df_prof['thickness_m'][0:idx+np.where(id_p)[0][0]])
                    Lrho = df_prof['density'][idx+np.where(id_p)[0][0]-1]
                    layfound = 2
                    if debug: print('[D]    Layfound value: ',layfound)
                else:
                    layfound = 0
                    if debug: print('[D]    Layfound value: ',layfound)
            else:
                layfound = 0
                
        if wlopt == 'nonper':
            
            id_p = selgt.str.contains('PP', regex=False) + selgt.str.contains('DF',regex = False) \
                        + selgt.str.contains('PPgp', regex=False)
            
            id_p = id_p.reset_index(drop=True)

            if sum(id_p.values)>0:
                if debug: print('[D]    WL grain type changed, looking for other non-persistent GT')
                Lupnew= df_prof[po][idx + np.where(id_p)[0][0]]
                if Lupnew > df_prof[po][-1:].values:
                    Llonew = df_prof[po][idx + np.where(id_p)[0][0]]
                    Lgtnew = df_prof['graintype'][idx + np.where(id_p)[0][0]]
                    SLdep = df_prof[po][0]-Lupnew        
                    SLrho = np.sum(df_prof['density'][0:(idx + np.where(id_p)[0][0])]  \
                                    * df_prof['thickness_m'][0:(idx + np.where(id_p)[0][0])]) \
                            / np.sum(df_prof['thickness_m'][0:(idx + np.where(id_p)[0][0])])
                    Lrho = df_prof['density'][(idx + np.where(id_p)[0][0])]
                    layfound = 2

                else:
                    layfound = 0 

            else:
                layfound = 0
                
    return Lupnew,Llonew,Lgtnew,Lrho,SLdep,SLrho,layfound


def initialize_with_RTA(config,season_list,df_met,df_P):
    """Initialize PAPs and DAPs with RTA of first (observed) profile
    -> Include date or index when to run this RTA search - could also be used for stations in between
    Input:
        season_list
        df_met
        df_P

    Output:
        df_P
    """

    """Thresholds"""
    debug       = int(config.get('AVAPRO', 'debug'))
    minSLdens   = int(config.get('AVAPRO-THOLDS-WL', 'minSLdens'))       # 120,min average density of a potential slab
    minSLthrsh  = float(config.get('AVAPRO-THOLDS-WL', 'minSLthrsh'))    # def.: .2,min. slab thickness [m] to be considered a healthy slab,

    ### Take first profile of list and iterate over its layers
    first_profile = season_list[0]
    # second_profile = season_list[1]

    RTA_last = 0
    RTA_thrsh = 0.7
    # pseudo_burial_days = 4
    per_gt = ['DH','SH','FC','FCxr']
    # nonper_gt = ['RG','MF','Mfcr']
    new_snow_gt = ['PP','DF', 'PPgp']
    pap_range_below_surf = 4 # check four layers below old snow surface 
    # burial_date = np.datetime64("NaT")
    burial_date = pd.Timestamp('1970-01-01')

    ### - Initialize either NAP or PAP - ###
    i_layer = 0
    while first_profile.loc[i_layer,'graintype'][0] in new_snow_gt:
        i_layer += 1
    
    if i_layer>0:
        # - snow new snow on older surface: check gt of old-surface #
        selgt_all = first_profile.graintype[i_layer:i_layer+pap_range_below_surf] # mabe implement chase where less than idx+3
        selgt = []
        for item in selgt_all:
            selgt.append(item[0])
        selgt = pd.Series(selgt)
        #% separate persistent and non-persistent
        id_p = selgt.str.contains('SH', regex=False)  + selgt.str.contains('FC',regex = False) \
            + selgt.str.contains('FCxr', regex=False) + selgt.str.contains('DH',regex =False)
        id_p = id_p.reset_index(drop=True)

        # - Only check old surfasce layer - #
        # gt_1 = first_profile.loc[i_layer,'graintype'][0] in per_gt
        # gt_2 = first_profile.loc[i_layer,'graintype'][1] in per_gt
        
        if sum(id_p.values) > 0:
            # - PAP found - #
            # get exact position of WL (uppermost persist. layer in window)
            po = 'height_m'
            df_P.loc[0, 'papburial'] = burial_date
            df_P.loc[0, 'papup']     = first_profile[po][((i_layer + id_p[id_p == True].index[0]))] #
            df_P.loc[0, 'paplo']     = first_profile[po][((i_layer + id_p[id_p == True].index[0])) + 1]
            df_P.loc[0, 'papgt'][:]  = first_profile.graintype[((i_layer + id_p[id_p == True].index[0]))]
            df_P.loc[0, 'papSLdep']  = first_profile[po][0] - df_P.loc[0, 'papup']

            rho = 'density'
            h  = 'thickness_m'
            df_P.loc[0, 'papSLrho']  = sum( first_profile[rho][0:i_layer  + id_p[id_p == True].index[0]] \
                                        * first_profile[h][0:i_layer  + id_p[id_p == True].index[0]])  \
                                        / sum(first_profile[h][0:i_layer  + id_p[id_p == True].index[0]])
            df_P.loc[0, 'papWLrho']  = first_profile[rho][i_layer  + id_p[id_p == True].index[0]]
            df_P.loc[0, 'papex']     = 1  # label persistent
            df_P.loc[0, 'napex']     = 0  # i.e. there is no non-persistent problem 
            if debug: print( '[D]   Initialization of pap: ', df_P.loc[0,'papgt'],'down',df_P.loc[0,'papSLdep']*100,'cm')
        else:
            # - NAP found - #
            rho = 'density'
            h = 'thickness_m' 
            df_P.loc[0, 'napburial'] = burial_date
            df_P.loc[0,'napup']      = first_profile['height_m'][i_layer]
            df_P.loc[0,'naplo']      = first_profile['height_m'][i_layer+1]
            df_P.loc[0,'napgt'][:]   = first_profile['graintype'][i_layer]
            df_P.loc[0,'napSLdep']   = first_profile['height_m'][0] - df_P['napup'][0]
            df_P.loc[0,'napSLrho']   =  sum( first_profile[rho][0:i_layer]* first_profile[h][0:i_layer]) / \
                                        sum(first_profile[h][0:i_layer])
            df_P.loc[0,'napWLrho']   = first_profile['density'][i_layer]
            df_P.loc[0,'napex']      = 1
            df_P.loc[0,'napcalc']    = 0

            # - need to evolve shear strength for the next day - #
            # tau_p, c_0 = fu_tau_p_CoJ15.fu_tau_p_CoJ15(0, df_P,  WLopt ='nap', opt = 'nonper') # %initial values defined there
            # df_P.loc[0,'napINI_tau_p'] = tau_p
            # df_P.loc[0,'napINI_c_0']   = c_0
            
            """Initialization of NAP"""
            index=0
            if df_P['napSLdep'][index] > minSLthrsh and df_P['napSLrho'][index] > minSLdens: # or df_P.drft[index] > drftthrsh: # Matlab 385
                if debug: print('[D]    Slab may ex. --> calc instab indices')
                df_P.loc[index,'napcalc'] = 1
                
                DAM,INI,DYN,PRO  = pro_instab.fu_instab(index, first_profile, df_P, df_met, WLopt='nap')
                df_P.loc[index,'napDAM_Sn']             = DAM[0] # natural stability index
                df_P.loc[index, 'napDAM_precstabMIN24'] = DAM[1]
                df_P.loc[index,'napDAM_extm2failMIN24'] = DAM[2]
                df_P.loc[index,'napDAM_tmcrit']         = DAM[3]                            
                    
                df_P.loc[index,'napINI_tau_p']          = INI[0]
                df_P.loc[index,'napINI_c_0']            = INI[1]
                df_P.loc[index,'napINI_Sk']             = INI[2]
                df_P.loc[index,'napINI_Sk_ana']         = INI[3] # % McClung & Sz 99, Monti etal 16
                df_P.loc[index,'napINI_mssANA']         = INI[4] 
                df_P.loc[index,'napINI_Sk_fem']         = INI[5] #% Rb & Sz 18 GRL append. 
                df_P.loc[index,'napINI_msswl']          = INI[6]
                df_P.loc[index,'napINI_sigma_g']        = INI[7]

                df_P.loc[index,'napDYN']                = DYN

                df_P.loc[index,'napPRO_ac_si06']        = PRO[0]                                 
                df_P.loc[index, 'napPRO_ac_vh16']       = PRO[1]                                 
                df_P.loc[index,'napPRO_wf']             = PRO[2]  
                df_P['napPRO_ac_ga17']                  = PRO[3]                     

            else:                     
                if debug: print('[D]    No healthy slab --> keep nap for +1 day')
                df_P.loc[index,'napcalc'] = 0
                #% need to evolve shear strength for the next day
                tau_p, c_0 = fu_tau_p_CoJ15.fu_tau_p_CoJ15(index, df_P,WLopt ='nap') # initial values defined there
                df_P.loc[index,'napINI_tau_p']  = tau_p
                df_P.loc[index,'napINI_c_0']    = c_0
            if debug: print('[D]  Initialization of NAP')

        i_old_surf = i_layer
    else:
        print('[I]  No NAP or PAP from initial snow profile')
        i_old_surf = 0
        pap_range_below_surf = 0

    """Initialization of DAPs"""
    if debug: print('[D]    Initialization of DAPs')
    DAP_i_layer = []
    DAP_up_ini  = []
    DAP_gt_ini  = []
    for i_layer in range(i_old_surf+pap_range_below_surf,len(first_profile)):
        gt_1 = first_profile.loc[i_layer,'graintype'][0] in per_gt
        gt_2 = first_profile.loc[i_layer,'graintype'][1] in per_gt
        if first_profile.loc[i_layer,'RTA'] >= RTA_thrsh and (gt_1 or gt_2):
            if first_profile.loc[i_layer,'RTA'] > RTA_last:
                # print(i_layer, first_profile.loc[i_layer,'height_m'], first_profile.loc[i_layer,'RTA'])
                DAP_up_ini.append(first_profile.loc[i_layer,'height_m'])
                DAP_gt_ini.append(first_profile.loc[i_layer,'graintype'])
                DAP_i_layer.append(i_layer)
                RTA_last = first_profile.loc[i_layer,'RTA'].copy()
        else:
            RTA_last = 0

    if len(DAP_up_ini) > 0:
        df_P.loc[0,'dapex'] = 1
    else:
        print('[I]  No DAP from initial snow profile')

    """Initialize all DAP parameters with values or NANs for further processing"""
    if df_P.loc[0, 'papex'] != 1:
        po = 'height_m'
        df_P.loc[0, 'papburial'] = burial_date
        df_P.loc[0, 'papup']    =  DAP_up_ini[0].copy()
        df_P.loc[0, 'papgt'][:] =  DAP_gt_ini[0].copy()
        df_P.loc[0, 'papex']    = 1
        df_P.loc[0, 'paplo']    =  first_profile[po][DAP_i_layer[0] + 1]
        df_P.loc[0, 'papSLdep'] = first_profile[po][0] - df_P.loc[0, 'papup']
        rho = 'density'
        h  = 'thickness_m'
        df_P.loc[0, 'papSLrho']  = sum( first_profile[rho][0: DAP_i_layer[0]] \
                                * first_profile[h][0:DAP_i_layer[0]])  \
                                / sum(first_profile[h][0:DAP_i_layer[0]])
        df_P.loc[0, 'papWLrho']  = first_profile[rho][DAP_i_layer[0]]
        if debug: print('[D]    Initialization of PAP from uppermost DAP: ', df_P.loc[0,'papgt'],'down',df_P.loc[0,'papSLdep']*100,'cm')

        df_P.loc[0,'dapup'][:] = DAP_up_ini[::-1][:-1].copy()
        #df_P.loc[0,'dapup'] = df_P.loc[0,'dapup'][0:-1].copy()
        df_P.loc[0,'dapgt'][:] = DAP_gt_ini[::-1][:-1].copy()
        #df_P.loc[0,'dapgt'][:] =df_P.loc[0,'dapgt'][0:-1].copy()
        df_P.loc[0,'dapburial'][:] = [burial_date] * (len(DAP_up_ini)-1)
        list_nan = [np.nan] * (len(DAP_up_ini)-1)
    else:
        df_P.loc[0,'dapup'][:] = DAP_up_ini[::-1].copy()
        df_P.loc[0,'dapgt'][:] = DAP_gt_ini[::-1].copy()
        df_P.loc[0,'dapburial'][:] = [burial_date] * len(DAP_up_ini)

        list_nan = [np.nan] * len(DAP_up_ini)

    dap_params = [  'daplo', 'dapSLdep', 'dapSLrho', 'dapWLrho',
                    'dapDAM_Sn','dapDAM_precstabMIN24','dapDAM_extm2failMIN24',
                    'dapDAM_tmcrit','dapINI_tau_p','dapINI_c_0','dapINI_Sk','dapINI_Sk_ana',
                    'dapINI_mssANA','dapINI_Sk_fem', 'dapINI_msswl','dapINI_sigma_g',
                    'dapDYN','dapPRO_ac_si06','dapPRO_ac_vh16','dapPRO_wf', 'dapPRO_ac_ga17']

    for param in dap_params:
        df_P.loc[0,param][:] = list_nan.copy()

    return df_P


def find_prior_problems(config,index,df_prof,df_met,df_P):
    """Move WLs from previous to current day, if still relevant
    Tasks:
        - update layer positions and properties of all DAPs via RTA initialization
    """

    """Thresholds"""
    debug       = int(config.get('AVAPRO', 'debug'))
    aggrgtDAPs  = int(config.get('AVAPRO-THOLDS-WL', 'aggrgtDAPs'))
    aggrgtdiff  = float(config.get('AVAPRO-THOLDS-WL', 'aggrgtdiff'))    # Aggregate DAPs if they are located within aggrgtdiff [m]
    outputDAP   = int(config.get('AVAPRO-THOLDS-WL', 'outputDAP'))       # Index of DAP to output; can create table: table(DAPLAYEROUT_txt,DAPLAYEROUT_dou)

    dropini     = int(config.get('AVAPRO-THOLDS-WL', 'dropini'))
    droppro     = float(config.get('AVAPRO-THOLDS-WL', 'droppro'))       # thresholds of initiation and propgation criteria for dropping (deep-)persitent weak layers
    drftthrsh   = float(config.get('AVAPRO-THOLDS-WL', 'drftthrsh'))     # min height of snow drifts [m], snp=:0.4m, cro:.2
    minSLdens   = int(config.get('AVAPRO-THOLDS-WL', 'minSLdens'))       # 120,min average density of a potential slab
    minnsthrsh  = float(config.get('AVAPRO-THOLDS-WL', 'minnsthrsh'))    # def.: 0.05m,min. amount of new snow [m], to start identifying problems (we find WL with windows, do don't care about +-5cm of fresh snow)
    nsthrsh     = float(config.get('AVAPRO-THOLDS-WL', 'nsthrsh'))       # def.: .2,amount of new snow [m] (within 48 or 24 hr), to keep a non-persis. problem for one more day
    minSLthrsh  = float(config.get('AVAPRO-THOLDS-WL', 'minSLthrsh'))    # def.: .2,min. slab thickness [m] to be considered a healthy slab,  

    vw_thrsh_days = int(config.get('AVAPRO-THOLDS-WL', 'vw_thrsh_days'))

    """Prior PAP"""
    if df_P.loc[index-1, 'papex'] == 1:
        if debug: print('[D]    Prior exiting PAP', index)   
        df_P.loc[index , 'papburial'] = df_P.loc[index-1 , 'papburial']
        #df_P['papup'][index], df_P['paplo'][index],df_P['papgt'][index], df_P['papWLrho'][index], df_P['papSLdep'][index],df_P['papSLrho'],layfound =
        Lupnew,Llonew,Lgtnew,Lrho,SLdep,SLrho,layfound = find_lay_prof(config,df_prof, df_P['papup'][index-1],df_P['papgt'][index-1], 'per')
        df_P.loc[index,'papup'] = Lupnew
        df_P.loc[index,'paplo'] = Llonew
        df_P.loc[index,'papgt'][:] = Lgtnew
        df_P.loc[index,'papWLrho'] = Lrho
        df_P.loc[index,'papSLdep'] = SLdep
        df_P.loc[index,'papSLrho'] = SLrho

        if layfound > 0: # matlab 152
        ### Treat problem, if persis. WL layer found
            df_P.loc[index,'papex']  = 1 # %keep the label "persistent"
            if (df_P['papSLdep'][index]> minSLthrsh and df_P['papSLrho'][index]>minSLdens) or df_P['drft'][index]>drftthrsh:
                if debug: print('[D]    Slab may exist, calculate instability for PAP', index)
                df_P.loc[index,'papcalc']  = 1

                DAM, INI, DYN, PRO = pro_instab.fu_instab(index,df_prof,df_P,df_met, WLopt='pap',pre_existing= True )     
                
                df_P.loc[index,'papDAM_Sn'] = DAM[0] # natural stability index
                df_P.loc[index,'papDAM_precstabMIN24'] = DAM[1]
                df_P.loc[index,'papDAM_extm2failMIN24']= DAM[2]
                df_P.loc[index, 'papDAM_tmcrit']  = DAM[3]
                df_P.loc[index, 'papINI_tau_p']  = INI[0]
                df_P.loc[index, 'papINI_c_0'] = INI[1]
                df_P.loc[index, 'papINI_Sk'] = INI[2]
                df_P.loc[index, 'papINI_Sk_ana'] = INI[3] # % McClung & Sz 99, Monti etal 16
                df_P.loc[ index, 'papINI_mssANA'] = INI[4] 
                df_P.loc[ index, 'papINI_Sk_fem'] = INI[5] #% Rb & Sz 18 GRL append. 
                df_P.loc[ index, 'papINI_msswl'] = INI[6]
                df_P.loc[ index, 'papINI_sigma_g'] = INI[7]
                df_P.loc[index, 'papDYN'] = DYN                    
                df_P.loc[index, 'papPRO_ac_si06']  = PRO[0]  
                df_P.loc[index, 'papPRO_ac_vh16']  = PRO[1]  
                df_P.loc[index, 'papPRO_wf']  = PRO[2]  
                df_P.loc[index, 'papPRO_ac_ga17']  = PRO[3]  

                if df_P.loc[index, 'papINI_tau_p'] /df_P.loc[index, 'papINI_mssANA'] > dropini \
                    or df_P.loc[index, 'papPRO_ac_vh16'] > droppro: # slcing
                    df_P.loc[index, 'papex']= 0
                    # print( ' INI:S=',(df_P['papINI'][index][0] / df_P['papINI'][index][4]*10)/10,' PROP:ac=',(df_P['papPRO'][index][1]*100),'cm / Sk_',alp,'=',(df_P['papINI'][index][2]*10)/10,' --> drop pap')
                #else: 
                #    print( ' INI:S=',(df_P['papINI_tau_p'][index] / df_P['papINI_mssANA'][index]*10)/10,' PROP:ac=',(df_P['papPRO_ac_vh16'][index]*100),'cm / Sk_',alp,'=',(df_P['papINI_Sk'][index]*10)/10)
            else:
                if debug: print('[D]    No slab, but keep PAP', index)
                df_P.loc[index, 'papcalc']  = 0 
                ### fu tau makes it slow; shear strength keeps evolving matlab 174
                tau_p, c_0 = fu_tau_p_CoJ15.fu_tau_p_CoJ15(index, df_P, WLopt ='pap', pre_existing= True) # rhowl only needed when layer born
                df_P.loc[index,'papINI_tau_p'] = tau_p
                df_P.loc[index,'papINI_c_0'] = c_0
        else:
            if debug: print('[D]    No WL found, drop PAP', index)
            df_P.loc[index,'papex']  = 2
            df_P.loc[index,'papcalc'] = 0           


    """Prior NAP"""
    if df_P.loc[index -1,'napex']  == 1:
        df_P.loc[index,'napburial'] = df_P.loc[index - 1,'napburial'] 
        Lupnew,Llonew,Lgtnew,Lrho,SLdep,SLrho,layfound = find_lay_prof(config,df_prof, df_P.loc[index-1, 'napup'],df_P.loc[index-1,'napgt'], 'nonper')
        df_P.loc[index,'napup']  = Lupnew
        df_P.loc[index,'naplo']  = Llonew
        df_P.loc[index,'napgt'][:]  = Lgtnew
        df_P.loc[index,'napWLrho']  = Lrho
        df_P.loc[index,'napSLdep']  = SLdep
        df_P.loc[index,'napSLrho']  = SLrho

        if layfound >0:

            ### Treat problem, if WL layer found
            if df_P.loc[index,'hn24' ] > nsthrsh or \
                (df_P.loc[index,'hn48' ] > 1* nsthrsh and df_P.loc[index,'hn24' ] > minnsthrsh ) or \
                (df_P.loc[index,'hn48']  > 1* nsthrsh and df_P.loc[index,'drft']  > drftthrsh):

                if debug: print('[D]    More snow added, investigate prior NAP again', index)
                df_P.loc[index, 'napcalc']= 1
                df_P.loc[index, 'napex'] = 1
                DAM, INI, DYN, PRO = pro_instab.fu_instab(index,df_prof,df_P,df_met, opt = 'nonper', WLopt='nap', pre_existing = True)     

                df_P.loc[index, 'napDAM_Sn']= DAM[0] # natural stability index
                df_P.loc[index,'napDAM_precstabMIN24'] = DAM[1]
                df_P.loc[index,'napDAM_extm2failMIN24'] = DAM[2]
                df_P.loc[index,'napDAM_tmcrit'] = DAM[3]                            
                
                df_P.loc[index,'napINI_tau_p']  = INI[0]
                df_P.loc[index,'napINI_c_0'] = INI[1]
                df_P.loc[index,'napINI_Sk'] = INI[2]
                df_P.loc[index,'napINI_Sk_ana'] = INI[3] # % McClung & Sz 99, Monti etal 16
                df_P.loc[index,'napINI_mssANA'] = INI[4] 
                df_P.loc[index,'napINI_Sk_fem'] = INI[5] #% Rb & Sz 18 GRL append. 
                df_P.loc[index,'napINI_msswl'] = INI[6]
                df_P.loc[index,'napINI_sigma_g'] = INI[7]
                
                df_P.loc[index,'napDYN'] = DYN

                df_P.loc[index,'napPRO_ac_si06'] = PRO[0]  
                df_P.loc[index,'napPRO_ac_vh16'] = PRO[1]  
                df_P.loc[index,'napPRO_wf'] = PRO[2]  
                df_P.loc[index,'napPRO_ac_ga17'] = PRO[3] 

            else:
                df_P.loc[index,'napcalc'] = 0
                df_P.loc[index,'napex'] = 0
        else:
            df_P.loc[index,'napex'] = 2
            df_P.loc[index,'napcalc'] = 0


    """Prior DAPs"""
    if df_P.loc[index - 1,'dapex'] == 1:
        ### DAP relevance (SLdep, SLrho): a dap could develop from an unrelevant pap

        ### Per def.: yesterday's problems are today's problems:
        ### Copy all prior deep problems to current day
        df_P.loc[index,'dapex'] = 1
        #df_P['dap'][index] = df_P['dap'][index-1] 

        df_P.loc[index,'dapup'][:] = df_P.loc[index-1,'dapup'].copy()
        df_P.loc[index,'daplo'][:] = df_P.loc[index-1,'daplo'].copy()
        df_P.loc[index,'dapgt'][:] = df_P.loc[index-1,'dapgt'].copy()
        df_P.loc[index,'dapWLrho'][:] = df_P.loc[index-1,'dapWLrho'].copy()
        df_P.loc[index,'dapSLdep'][:] = df_P.loc[index-1,'dapSLdep'].copy()
        df_P.loc[index,'dapSLrho'][:] = df_P.loc[index-1,'dapSLrho'].copy()
        df_P.loc[index,'dapburial'][:] = df_P.loc[index-1,'dapburial'].copy()
                        
        df_P.loc[index,'dapDAM_Sn'][:] = df_P.loc[index-1,'dapDAM_Sn'].copy()
        df_P.loc[index,'dapDAM_precstabMIN24'][:] = df_P.loc[index-1,'dapDAM_precstabMIN24'].copy()
        df_P.loc[index,'dapDAM_extm2failMIN24'][:] = df_P.loc[index-1,'dapDAM_extm2failMIN24'].copy() 
        df_P.loc[index,'dapDAM_tmcrit'][:] = df_P.loc[index-1,'dapDAM_tmcrit'].copy()
                        
        df_P.loc[index,'dapINI_tau_p'][:] = df_P.loc[index-1,'dapINI_tau_p'].copy()
        df_P.loc[index,'dapINI_c_0'][:] = df_P.loc[index-1,'dapINI_c_0'].copy()
        df_P.loc[index,'dapINI_Sk'][:] = df_P.loc[index-1,'dapINI_Sk'].copy()
        df_P.loc[index,'dapINI_Sk_ana'][:] = df_P.loc[index-1,'dapDAM_Sn'].copy()
        df_P.loc[index,'dapINI_mssANA'][:] = df_P.loc[index-1,'dapINI_mssANA'].copy()
        df_P.loc[index,'dapINI_Sk_fem'][:] = df_P.loc[index-1,'dapINI_Sk_fem'].copy()              
        df_P.loc[index,'dapINI_msswl'][:] = df_P.loc[index-1,'dapINI_msswl'].copy()
        df_P.loc[index,'dapINI_sigma_g'][:] = df_P.loc[index-1,'dapINI_sigma_g'].copy()
                    
        df_P.loc[index,'dapPRO_ac_si06'][:]  = df_P.loc[index-1,'dapPRO_ac_si06'].copy()
        df_P.loc[index,'dapPRO_ac_vh16'][:] = df_P.loc[index-1,'dapPRO_ac_vh16'].copy()
        df_P.loc[index,'dapPRO_wf'][:] = df_P.loc[index-1,'dapPRO_wf'].copy()
        df_P.loc[index,'dapPRO_ac_ga17'][:]  = df_P.loc[index-1,'dapPRO_ac_ga17'].copy()

        df_P.loc[index,'dapDYN'][:]  = df_P.loc[index-1,'dapDYN'].copy()
        
        for ele in range(len(df_P['dapup'][index])):
            ### No processing of dropped daps or daps settling towards the ground
            if ~np.isnan( df_P.loc[index-1, 'dapup'][ele]) and  df_P.loc[index-1, 'dapup'][ele] > 0.1 :
                #if df_P['dapgt'][index-1][ele][0] == 'FCxr':  ## for matlb comparison (matlab loses FCxr)
                    #   df_P['dapgt'][index-1][ele][0] = 'FC'
                  
                Lupnew,Llonew,Lgtnew,Lrho,SLdep,SLrho,layfound = find_lay_prof(config,df_prof, df_P.loc[index-1,'dapup'][ele],df_P.loc[index-1,'dapgt'][ele], 'per') 
                df_P.loc[index,'dapup'][ele] = Lupnew
                df_P.loc[index,'daplo'][ele] = Llonew
                
                df_P.loc[index,'dapWLrho'][ele] = Lrho
                df_P.loc[index,'dapSLdep'][ele] = SLdep
                df_P.loc[index,'dapSLrho'][ele] = SLrho

                if layfound > 0:
                    DAM, INI, DYN, PRO = pro_instab.fu_instab(index,df_prof,df_P,df_met, WLopt='dap', dap_nr = ele)   
                    
                    df_P.loc[index,'dapDAM_Sn'][ele] = DAM[0]
                    df_P.loc[index,'dapDAM_precstabMIN24'][ele] = DAM[1]
                    df_P.loc[index,'dapDAM_extm2failMIN24'][ele] = DAM[2]
                    df_P.loc[index,'dapDAM_tmcrit'][ele] = DAM[3]
                        
                    df_P.loc[index,'dapINI_tau_p'][ele] = INI[0]
                    df_P.loc[index,'dapINI_c_0'][ele] = INI[1]
                    df_P.loc[index,'dapINI_Sk'][ele] = INI[2]
                    df_P.loc[index,'dapINI_Sk_ana'][ele] = INI[3]
                    df_P.loc[index,'dapINI_mssANA'][ele] = INI[4]
                    df_P.loc[index,'dapINI_Sk_fem'][ele] = INI[5]                              
                    df_P.loc[index,'dapINI_msswl'][ele] = INI[6]
                    df_P.loc[index,'dapINI_sigma_g'][ele] = INI[7]
                    
                    df_P.loc[index,'dapPRO_ac_si06'][ele] = PRO[0]
                    df_P.loc[index,'dapPRO_ac_vh16'][ele] = PRO[1]
                    df_P.loc[index,'dapPRO_wf'][ele] = PRO[2]
                    df_P.loc[index,'dapPRO_ac_ga17'][ele] = PRO[3]

                    df_P.loc[index,'dapDYN'][ele] = DYN

                    if df_P.loc[index,'dapINI_tau_p'][ele]/df_P.loc[index,'dapINI_mssANA'][ele]> dropini or df_P.loc[index,'dapPRO_ac_vh16'][ele]> droppro:
                        df_P.loc[index,'dapup'][ele] = np.NaN #  % that means dropping, see ifloop
                        #disp(['    INI:S=',num2str(round(P.dap(mm).INI(nn,1)/P.dap(mm).INI(nn,5)*10)/10),' PROP:ac=',num2str(round(P.dap(mm).PRO(nn,2)*100)),'cm / Sk_',num2str(alp),'=',num2str(round(P.dap(mm).INI(nn,3)*10)/10),' --> drop dap'])
                        if debug: print('[D]    Drop DAP', index)
                    else:
                        #disp(['    INI:S=',num2str(round(P.dap(mm).INI(nn,1)/P.dap(mm).INI(nn,5)*10)/10),' PROP:ac=',num2str(round(P.dap(mm).PRO(nn,2)*100)),'cm / Sk_',num2str(alp),'=',num2str(round(P.dap(mm).INI(nn,3)*10)/10)])
                        if debug: print('[D]    DAP still relevant', index)
                else:
                    if debug: print('[D]    Drop DAP (No WL found)', index)
                    df_P.loc[index,'dapup'][ele] = np.NaN
                    
        if aggrgtDAPs == 1 and sum( ~np.isnan(df_P['dapup'][index]))>2: 
            ind_arg = np.NaN
            for i,x in enumerate(np.diff(df_P['dapup'][index][::])) :
                if x < aggrgtdiff:
                    ind_arg = i
                    
            if ~np.isnan(ind_arg):
                # ind_arg = len(df_P['dapINI_tau_p'][index])-ind_arg-1
                # print (index, df_P['dapINI_tau_p'][index])
                # print(' ')
                # print(ind_arg)
                b = np.argmax(df_P['dapINI_tau_p'][index][ind_arg:ind_arg+2])
                ind_arg = ind_arg + b      
                df_P.loc[index,'dapup'][ind_arg] = np.NaN #% that means dropping, see ifloop
                if debug: print('[D]    Drop DAP', index)
        if outputDAP !=0:
            print('[D]    Implement matlab line 254')


    """WSAPs (Wind slab avalanche problem)"""
    if (df_P.loc[index-1, 'winex'] == 1) and (df_P.loc[index-1, 'winex_count'] < vw_thrsh_days):
        df_P.loc[index, 'winex'] = 1
        df_P.loc[index, 'winex_count'] = df_P.loc[index-1, 'winex_count'] + 1
    else:
        df_P.loc[index, 'winex_count'] = 0
    
    return df_P


def identify_wet_ap(config,index,ind_wet,season_list_red,df_met_red,df_P,cycle_nr):
    """WAP"""

    """Thresholds"""
    debug       = int(config.get('AVAPRO', 'debug'))
    minSLthrsh  = float(config.get('AVAPRO-THOLDS-WL', 'minSLthrsh'))    # def.: .2,min. slab thickness [m] to be considered a healthy slab,
    lwcthrsh_0  = float(config.get('AVAPRO-THOLDS-WL', 'lwcthrsh_0'))    # Christophs LWC threshold for maztrix flow: 0.03; threshold for preferential flow: 0.01
    lwcthrsh_1  = float(config.get('AVAPRO-THOLDS-WL', 'lwcthrsh_1'))

    # % proxy for crit water: [kg] total water in snow cover at wettime > gutdünken * mass of snow cover?
    df_P.loc[index ,'wapTWAT'] = df_met_red.loc[ind_wet, 'MS_Water']
    df_P.loc[index,'wapSWE'] = df_met_red.loc[ind_wet , 'SWE']
    
    if df_met_red.MS_Water.iloc[ind_wet] > 0.0002 * df_met_red.SWE.iloc[ind_wet]: # %gutdünken=20%?? (für snp_WFJ201718:2%??)
    # if df_met_red.MS_Water.iloc[ind_wet] >= 0 * df_met_red.SWE.iloc[ind_wet]:
        #sel timestamp of wettime
        dy = df_met_red.loc[ind_wet , 'timestamp']
        
        ### Get profile from wet time
        df_prof_wet = season_list_red[ind_wet] 
        if df_prof_wet.date[0] != dy:
            if debug: print('[D]    Date of season_list_red (profiles) -',df_prof_wet.date[0],'- does not match date of df_P (wet cycle)-',dy)

        ### Check if all layers are wet except ground (---> isothermal) line 279
        df_prof_wet['lwc_tot'] = df_prof_wet['lwc'] /100
        df_P.loc[index, 'wapISO'] = df_prof_wet[df_prof_wet['lwc_tot']>0].shape[0] >= (df_prof_wet['lwc_tot'].shape[0]-1)
    
        ### Ex. still a dry snow layer in snow cover? (except ground)
        po= 'height_m'# just for better readability
        wet_po = df_prof_wet[po][:-1] < df_prof_wet[po][0] - minSLthrsh
        wet_lwc = df_prof_wet['lwc_tot'][:-1]<0.011
        gt_1 = df_prof_wet.graintype[:-1].str[0].str.contains('MF', case = False )
        # gt_2 = df_prof_wet.graintype[:-1].str[1].str.contains( 'MF', case = False ) not needed at the moment
        df_P.loc[index, 'wapWLdry'] = sum( wet_po & wet_lwc & gt_1)>0 #ex. still a dry snow layer in snow cover? (except ground)

        ### [1] liquid water content as volume fraction lwctdysmet = M.totwater; % compare with smet file -- total liquid wate
        h = 'thickness_m'
        lwctdy =sum( df_prof_wet['lwc_tot']*df_prof_wet[h]) / sum(df_prof_wet[h])

        ### Write LWC to snow instab. data 
        df_P.loc[index, 'wapLWC'] = lwctdy     
        # print('  wet now? ... ', df_met_red.MS_Water.iloc[indwet]/df_met_red.SWE.iloc[indwet],'[m]','lwc:', lwctdy,'[v]') # einheiten kontrolieren
        
        ### Track of wetcycle and assigne threshold according to Mitterer et al 2016
        if cycle_nr >= 2 and df_P.loc[index-1, 'wapcycle'] != 1:
            lwc_thrsh = lwcthrsh_1
        else: ### cycle_nr == 1
            lwc_thrsh = lwcthrsh_0

        if lwctdy > lwc_thrsh:
            # print('-->wap')
            df_P.loc[index,'wapex']= 1 # % wet snow problem ex.
            if ~np.isnan(df_P.loc[index -1 , 'wapex']):
                df_P.loc[index, 'waponset'] = df_P.loc[index-1, 'waponset']
                df_P.loc[index, 'wapcycle'] = df_P.loc[index-1, 'wapcycle']
            else:
                df_P.loc[index, 'waponset'] = dy
                df_P.loc[index, 'wapcycle'] = cycle_nr
                cycle_nr =  cycle_nr + 1
        else:
            if debug: print('[D]    No WAP') 
    else:
        df_P.loc[index, 'wapex']= np.NaN
    
    return df_P, cycle_nr


def identify_new_nap_or_pap(config,index,ind_dry,season_list_red,df_met_red,df_met,df_P,deltahs):
    """Identify new WL at old surface (surface of previous profile) as NAP or PAP"""

    """Thresholds"""
    debug       = int(config.get('AVAPRO', 'debug'))
    drftthrsh   = float(config.get('AVAPRO-THOLDS-WL', 'drftthrsh'))     # min height of snow drifts [m], snp=:0.4m, cro:.2
    minSLdens   = int(config.get('AVAPRO-THOLDS-WL', 'minSLdens'))       # 120,min average density of a potential slab
    minSLthrsh  = float(config.get('AVAPRO-THOLDS-WL', 'minSLthrsh'))    # def.: .2,min. slab thickness [m] to be considered a healthy slab,
    minnsthrsh  = float(config.get('AVAPRO-THOLDS-WL', 'minnsthrsh'))    # def.: 0.05m,min. amount of new snow [m], to start identifying problems (we find WL with windows, do don't care about +-5cm of fresh snow)

    vw_thrsh = int(config.get('AVAPRO-THOLDS-WL', 'vw_thrsh'))

    po = 'height_m'
    if index == 0:
        pass
    elif df_P.loc[index-1, 'napex'] !=1  and df_P.loc[index, 'hn24'] > minnsthrsh: # matlab 304
        dy = df_met_red.loc[ind_dry , 'timestamp']
        df_prof = season_list_red[ind_dry]
        if df_prof.date[0] != dy:
            print('napex, date of season_list_red (profiles) -',df_prof.date[0],'- does not match date of df_P (wet cycle)-',dy)
        
        idx = sum((df_prof[po]-(df_prof[po][0] - deltahs))>0) #%index to old surface !! hn24 used instead before
        if idx > 0 and idx < (len(df_prof[po])-2): 
            if debug: print("[D]    Potential new WL: ",index, df_P.dy[index])
            ### Window to check out properties of old snow surface, get grain types "around" old surface, min fun prevents so exceed length of array
            selgt_all = df_prof.graintype[idx-1:idx+4]
            selgt = []
            for item in selgt_all:
                selgt.append(item[0])
            selgt = pd.Series(selgt)
            
            ### Separate persistent and non-persistent
            id_p = selgt.str.contains('SH', regex=False)  + selgt.str.contains('FC',regex = False) \
                + selgt.str.contains('FCxr', regex=False) + selgt.str.contains('DH',regex =False)
            id_p = id_p.reset_index(drop=True)

            if sum(id_p.values)> 0:  # matlab 319       
                """New PAP, old PAP is moved to DAPs"""
                if df_P.loc[index, 'papex'] == 1:
                    if debug: print("[D]    New PAP, old PAP becomes DAP: ",df_P.loc[index, 'dy'])
                    df_P.loc[index, 'dapex'] = 1 

                    df_P.loc[index,'dapburial'].append(str(df_P['papburial'][index]))
                    df_P.loc[index,'dapup'].append(df_P['papup'][index])                                             
                    df_P.loc[index, 'daplo'].append(df_P['paplo'][index])
                    
                    df_P.loc[index,'dapgt'].append(df_P['papgt'][index]) 
                    df_P.loc[index,'dapSLdep'].append(df_P['papSLdep'][index])
                    df_P.loc[index,'dapSLrho'].append(df_P['papSLrho'][index])
                    df_P.loc[index,'dapWLrho'].append(df_P['papWLrho'][index])
                    
                    df_P.loc[index,'dapDAM_Sn'].append(df_P['papDAM_Sn'][index] ) 
                    df_P.loc[index,'dapDAM_precstabMIN24'].append( df_P['papDAM_precstabMIN24'][index])
                    df_P.loc[index,'dapDAM_extm2failMIN24'].append( df_P['papDAM_extm2failMIN24'][index])
                    df_P.loc[index,'dapDAM_tmcrit'].append(df_P['papDAM_tmcrit'][index])
                                                        
                    df_P.loc[index,'dapINI_tau_p'].append(df_P['papINI_tau_p'][index])
                    df_P.loc[index,'dapINI_c_0'].append(df_P['papINI_c_0'][index])
                    df_P.loc[index,'dapINI_Sk'].append(df_P.loc[index,'papINI_Sk'])
                    df_P.loc[index,'dapINI_Sk_ana'].append(df_P.loc[index,'papINI_Sk_ana'])
                    df_P.loc[index,'dapINI_mssANA'].append(df_P.loc[index,'papINI_mssANA'])
                    df_P.loc[index,'dapINI_Sk_fem'].append(df_P.loc[index,'papINI_Sk_fem'])                                
                    df_P.loc[index,'dapINI_msswl'].append(df_P.loc[index,'papINI_msswl'])
                    df_P.loc[index,'dapINI_sigma_g'].append(df_P.loc[index,'papINI_sigma_g'])
                                                    
                    df_P.loc[index,'dapPRO_ac_si06'].append(df_P.loc[index,'papPRO_ac_si06'])
                    df_P.loc[index,'dapPRO_ac_vh16'].append(df_P.loc[index,'papPRO_ac_vh16'])
                    df_P.loc[index,'dapPRO_wf'].append(df_P.loc[index,'papPRO_wf'])
                    df_P.loc[index,'dapPRO_ac_ga17'].append(df_P.loc[index,'papPRO_ac_ga17'])
                    
                    df_P.loc[index,'dapDYN'].append(df_P.loc[index,'papDYN'])
                    
                ### reset the pap instability properties (have just been written to dap)
                df_P.loc[index,'papDAM_Sn']= np.NaN # natural stability index
                df_P.loc[index,'papDAM_precstabMIN24'] = np.NaN
                df_P.loc[index,'papDAM_extm2failMIN24'] = np.NaN
                df_P.loc[index,'papDAM_tmcrit'] = np.NaN

                df_P.loc[index,'papINI_tau_p'] = np.NaN
                df_P.loc[index,'papINI_c_0']  = np.NaN
                df_P.loc[index,'papINI_Sk']  = np.NaN
                df_P.loc[index,'papINI_Sk_ana']  = np.NaN # % McClung & Sz s99, Monti etal 16
                df_P.loc[index,'papINI_mssANA']  = np.NaN 
                df_P.loc[index,'papINI_Sk_fem']  = np.NaN #% Rb & Sz 18 GRL append. 
                df_P.loc[index,'papINI_msswl']  = np.NaN
                df_P.loc[index,'papINI_sigma_g']  = np.NaN

                df_P.loc[index,'papDYN']  = np.NaN
                df_P.loc[index,'papPRO_ac_si06']  = np.NaN  
                df_P.loc[index,'papPRO_ac_vh16'] = np.NaN  
                df_P.loc[index,'papPRO_wf']  = np.NaN  
                df_P.loc[index,'papPRO_ac_ga17']  = np.NaN

                # get exact position of WL (uppermost persist. layer in window)
                df_P.loc[index, 'papburial'] = df_met_red.timestamp.iloc[ind_dry] # ???
                df_P.loc[index, 'papup'] = df_prof[po][((idx-1 + id_p[id_p == True].index[0]))] #
                df_P.loc[index, 'paplo'] = df_prof[po][((idx-1 + id_p[id_p == True].index[0])) + 1]
                df_P.loc[index, 'papgt'][:] =  df_prof.graintype[((idx-1 + id_p[id_p == True].index[0]))]
                df_P.loc[index, 'papSLdep'] = df_prof[po][0] - df_P.papup[index]

                rho = 'density'
                h  = 'thickness_m'
                
                df_P.loc[index, 'papSLrho'] =   sum( df_prof[rho][0:idx  + id_p[id_p == True].index[0]] \
                                                * df_prof[h][0:idx  + id_p[id_p == True].index[0]])  \
                                                / sum(df_prof[h][0:idx  + id_p[id_p == True].index[0]])
                df_P.loc[index, 'papWLrho'] = df_prof[rho][idx + id_p[id_p == True].index[0]]
                df_P.loc[index, 'papex']    = 1  # label persistent
                df_P.loc[index, 'napex']    = 0  # i.e. there is no non-persistent problem 
                if debug: print('[D]    New PAP: ', df_P.papgt[index],'down',[df_P.papSLdep[index]*100],'cm')

                ### Treat problem (Matlab 353)
                if (df_P.loc[index, 'papSLdep'] > minSLthrsh and df_P.loc[index, 'papSLrho'] > minSLdens) \
                    or df_P.loc[index, 'drft'] > drftthrsh :
                    #  slab may ex. --> calculate INI and PROP
                    df_P.loc[index, 'papcalc'] = 1
                    DAM, INI, DYN, PRO = pro_instab.fu_instab(index,df_prof,df_P,df_met, WLopt='pap')     

                    df_P.loc[index, 'papDAM_Sn'] = DAM[0] # natural stability index
                    df_P.loc[index,'papDAM_precstabMIN24'] = DAM[1]
                    df_P.loc[index,'papDAM_extm2failMIN24']  = DAM[2]
                    df_P.loc[index,'papDAM_tmcrit'] = DAM[3]

                    df_P.loc[index,'papINI_tau_p'] = INI[0]
                    df_P.loc[index,'papINI_c_0']  = INI[1]
                    df_P.loc[index,'papINI_Sk']  = INI[2]
                    df_P.loc[index,'papINI_Sk_ana']  = INI[3] # % McClung & Sz 99, Monti etal 16
                    df_P.loc[index,'papINI_mssANA']  = INI[4] 
                    df_P.loc[index,'papINI_Sk_fem']  = INI[5] #% Rb & Sz 18 GRL append. 
                    df_P.loc[index,'papINI_msswl']  = INI[6]
                    df_P.loc[index,'papINI_sigma_g']  = INI[7]

                    df_P.loc[index, 'papDYN']  = DYN

                    df_P.loc[index,'papPRO_ac_si06'] = PRO[0]  
                    df_P.loc[index,'papPRO_ac_vh16'] = PRO[1]  
                    df_P.loc[index,'papPRO_wf'] = PRO[2]  
                    df_P.loc[index,'papPRO_ac_ga17'] = PRO[3]  

                else : ### matlab 363
                    if debug: print('[D]    No slab yet, but PAP kept')
                    df_P.loc[index, 'papcalc'] = 0
                    tau_p, c_0 = fu_tau_p_CoJ15.fu_tau_p_CoJ15(index, df_P) 
                    df_P.loc[index,'papINI_tau_p'] = tau_p
                    df_P.loc[index,'papINI_c_0']  = c_0

            else: ### matlab 371
                """New NAP, its position is the old snow surface"""
                df_P.loc[index, 'napburial'] = df_met_red['timestamp'][ind_dry]
                df_P.loc[index,'napup'] = df_prof['height_m'][idx]
                df_P.loc[index,'naplo'] = df_prof['height_m'][idx+1]
                df_P.loc[index,'napgt'][:] = df_prof['graintype'][idx]
                df_P.loc[index,'napSLdep'] = df_prof['height_m'][0] - df_P['napup'][index]
                rho = 'density'
                h = 'thickness_m' 

                df_P.loc[index,'napSLrho'] =  sum( df_prof[rho][0:idx+1]* df_prof[h][0:idx+1]) / \
                                              sum(df_prof[h][0:idx+1])
                df_P.loc[index,'napWLrho'] = df_prof['density'][idx]
                df_P.loc[index,'napex'] = 1 # %label non-persistent
                
                ### Treat or drop problem  
                if df_P['napSLdep'][index] > minSLthrsh and df_P['napSLrho'][index] > minSLdens \
                    or df_P.drft[index] > drftthrsh: # Matlab 385

                    df_P.loc[index,'napcalc'] = 1
                    DAM,INI,DYN,PRO  = pro_instab.fu_instab(index, df_prof, df_P, df_met, WLopt='nap')
                    df_P.loc[index,'napDAM_Sn'] = DAM[0] # natural stability index
                    df_P.loc[index,'napDAM_precstabMIN24'] = DAM[1]
                    df_P.loc[index,'napDAM_extm2failMIN24'] = DAM[2]
                    df_P.loc[index,'napDAM_tmcrit'] = DAM[3]                            
                        
                    df_P.loc[index,'napINI_tau_p'] = INI[0]
                    df_P.loc[index,'napINI_c_0'] = INI[1]
                    df_P.loc[index,'napINI_Sk'] = INI[2]
                    df_P.loc[index,'napINI_Sk_ana'] = INI[3] # % McClung & Sz 99, Monti etal 16
                    df_P.loc[index,'napINI_mssANA'] = INI[4] 
                    df_P.loc[index,'napINI_Sk_fem'] = INI[5] #% Rb & Sz 18 GRL append. 
                    df_P.loc[index,'napINI_msswl'] = INI[6]
                    df_P.loc[index,'napINI_sigma_g'] = INI[7]

                    df_P.loc[index,'napDYN'] = DYN
                
                    df_P.loc[index,'napPRO_ac_si06'] = PRO[0]                                 
                    df_P.loc[index, 'napPRO_ac_vh16'] = PRO[1]                                 
                    df_P.loc[index,'napPRO_wf'] = PRO[2]  
                    df_P['napPRO_ac_ga17'] = PRO[3]                     

                else:   # Matlab line 392                     
                    if debug: print('[D]    No healthy slab, keep NAP for +1 day')
                    df_P.loc[index,'napcalc'] = 0
                    tau_p, c_0 = fu_tau_p_CoJ15.fu_tau_p_CoJ15(index, df_P, WLopt ='nap')
                    df_P.loc[index,'napINI_tau_p']  = tau_p
                    df_P.loc[index,'napINI_c_0']  = c_0
        
        """Wind only case"""
        ### Matlab line 401
        ### Always considered. Also when surface does not change     
        if df_P.napcalc[index] != 1 and df_P.drft.iloc[index] > drftthrsh:
            print("!!!! Wind only entered !!!!")
            ### In this case flatfield sim do not make a slab --> make something up: see fu_instab opt=windonly
            df_P.loc[index, 'napwindonly'] = 1
            df_P.loc[index, 'napex'] = 2 #  % no weak layer -> we can't pull up the problem the following day
            DAM, INI, DYN, PRO = pro_instab.fu_instab(index, df_prof, df_P,df_met, opt ='windonly', WLopt='nap', winddrift= True ) # initial values defined there
            
            df_P.loc[index,'napDAM_Sn'] = DAM[0] # natural stability index
            df_P.loc[index,'napDAM_precstabMIN24']  = DAM[1]
            df_P.loc[index,'napDAM_extm2failMIN24'] = DAM[2]
            df_P.loc[index,'napDAM_tmcrit'] = DAM[3]
            
            df_P.loc[index,'napINI_tau_p'] = INI[0]
            df_P.loc[index,'napINI_c_0'] = INI[1]
            df_P.loc[index,'napINI_Sk'] = INI[2]
            df_P.loc[index,'napINI_Sk_ana'] = INI[3] # McClung & Sz 99, Monti etal 16
            df_P.loc[index,'napINI_mssANA'] = INI[4] 
            df_P.loc[index,'napINI_Sk_fem'] = INI[5] # Rb & Sz 18 GRL append. 
            df_P.loc[index,'napINI_msswl'] = INI[6]
            df_P.loc[index,'napINI_sigma_g'] = INI[7]
            df_P.loc[index,'napDYN'] = DYN
            
            #if not pd.isna(PRO.all()):  
            df_P.loc[index,'napPRO_ac_si06']  = PRO[0]  
            df_P.loc[index,'napPRO_ac_vh16'] = PRO[1]  
            df_P.loc[index,'napPRO_wf']  = PRO[2]  
            df_P.loc[index,'napPRO_ac_ga17']  = PRO[3] 
            print('[D]  Wind slabs possibly /wo new snow')

        """Wind only based on Meteo df and surface snow (cmb)"""
        ### Wind slab problem if wind above threshold and loose snow on top within last (3) days
        ### index: index in df_P
        ### ind_dry and ind_wet: index in season_list_red
        dy_past     = dy-datetime.timedelta(days=1)
        mask_wind   = (df_met['timestamp'] >= dy_past) & (df_met['timestamp'] <= dy)
        df_met_wind = df_met.loc[mask_wind].reset_index(drop=True)
        mask_vw_threshold = df_met_wind['VW_drift'] >= vw_thrsh
        # print(df_met_wind["timestamp"][mask_vw_threshold])
        
        if np.any(mask_vw_threshold):
            df_prof_wind = season_list_red[ind_dry]
        
            ### Check surface grain type (if not 'RG' snow can be transported)
            if df_prof_wind['graintype'][0][0] in ['PP','PPgp','DF','FC','SH','DH']:
                df_P.loc[index, 'winex'] = 1
                df_P.loc[index, 'winex_count'] = 1

    return df_P