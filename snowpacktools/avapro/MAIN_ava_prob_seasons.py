##### Matlab to python project
# ### imports ###
import os
import sys
import pandas as pd
import numpy as np
import pickle
import glob
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.ticker import AutoMinorLocator, FuncFormatter
import datetime
import configparser


#import processor_PRO as pro_pro  # function to read in .pro files
# from snowpro import snowpro 
# import read_smet_file as pro_smet # functions to read in meteo data from .smet file
import find_aps as pro_find_aps # skript to detect weaklayers and stability calc
import post_processing as pro_post_processing # assign avaprobs to relevant wl
#config_file = '/Users/martin/Documents/LWD_extension/avaproblemsmatlab2python/py-code/avaprob.ini'

def main(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)
    print('...reading in config parameters from ini file...')
    # get litst of input files
    rerun_find_WL = int(config.get('avaProbs_general', 'rerun_find_WL')) 
    rerun_assign_avaprobs = int(config.get('avaProbs_general', 'rerun_assign_avaprobs')) 
    initilization_type = config.get('avaProbs_general', 'initilization_type') 
    SNP_INPUT_FOLDER = config.get('avaProbs_general', 'SNP_INPUT_FOLDER')

    list_pro = sorted(glob.glob(SNP_INPUT_FOLDER + "/*.pro"))
    list_smet = sorted(glob.glob(SNP_INPUT_FOLDER + "/*.smet"))
    list_csv =  sorted(glob.glob(SNP_INPUT_FOLDER + "/*.csv"))
    print(list_pro, 'available in INPUT folder...')
    # Reading in the parameters from ini file
        # -- PARAMETER & METADATA -- # definend in ini file
    station_name = SNP_INPUT_FOLDER.split('/')[-2]
    output_path = config.get('avaProbs_general', 'output_path') 
    path = output_path + station_name

    # getting aspects for simulation
    slope_aspects = config.get('avaProbs_general', 'SLOPE_ASPECTS').split(' ')

    slope_aspect_dict= {'flat' : station_name,
                        'N': station_name+'1',
                        'E' : station_name+'2',
                        'S' : station_name+'3',
                        'W' : station_name+'4',
                        }
    slope_aspects_sim = []
    for aspect in slope_aspects:
        slope_aspects_sim.append(slope_aspect_dict[aspect])
    print(slope_aspects, '...aspects for simulations selected...')
    #reduce pro and smet list to aspects to simulate:
    list_pro_red = []
    list_smet_red = []
    for aspects in slope_aspects_sim:
        for ele in list_pro:
            element = ele.split('/')
            element = element[-1].split('.')[0]
            if aspects == element:
                list_pro_red.append(ele)
                print(ele, 'added to simulation')
                break
        for ele in list_smet:
            element = ele.split('/')
            element = element[-1].split('.')[0]
            if aspects == element:
                list_smet_red.append(ele)
                print(ele, 'added to simulation')
                break
    #check if .pro has same len as .smet:
    if len(list_pro_red) == len(list_smet_red):
        pass
    elif len(list_pro_red) == 0:
        raise ValueError('no input, please check...')
    else:
        raise ValueError('.pro input does not have the same len as .smet. check your input')

    # defining season
    season_start = config.get('avaProbs_general', 'season_start') 
    season_end = config.get('avaProbs_general', 'season_end') 
    drytime=int(config.get('avaProbs_general', 'drytime')) 
    wettime=int(config.get('avaProbs_general', 'wettime')) 
    scmopt = config.get('avaProbs_general', 'scmopt')
    release = config.get('avaProbs_general', 'release')

    #thresholds for detection of WL
    thresholds_aps = {}
    thresholds_aps['alp'] = int(config.get('avaProbs_thresholds_aps', 'alp'))
    thresholds_aps['calcFEM'] = int(config.get('avaProbs_thresholds_aps', 'calcFEM'))
    thresholds_aps['owSCtaup'] = int(config.get('avaProbs_thresholds_aps', 'owSCtaup'))
    thresholds_aps['lookat'] = int(config.get('avaProbs_thresholds_aps', 'lookat'))
    thresholds_aps['aggrgtDAPs'] = int(config.get('avaProbs_thresholds_aps', 'aggrgtDAPs'))
    thresholds_aps['aggrgtdiff'] = float(config.get('avaProbs_thresholds_aps', 'aggrgtdiff'))
    thresholds_aps['outputDAP'] = int(config.get('avaProbs_thresholds_aps', 'outputDAP'))
    thresholds_aps['dropini'] = int(config.get('avaProbs_thresholds_aps', 'dropini'))
    thresholds_aps['droppro'] = float(config.get('avaProbs_thresholds_aps', 'droppro'))
    thresholds_aps['drftthrsh'] = float(config.get('avaProbs_thresholds_aps', 'drftthrsh'))
    thresholds_aps['minSLdens'] = int(config.get('avaProbs_thresholds_aps', 'minSLdens'))
    thresholds_aps['minnsthrsh'] = float(config.get('avaProbs_thresholds_aps', 'minnsthrsh'))
    thresholds_aps['nsthrsh'] = float(config.get('avaProbs_thresholds_aps', 'nsthrsh'))
    thresholds_aps['minSLthrsh'] = float(config.get('avaProbs_thresholds_aps', 'minSLthrsh'))
    thresholds_aps['lwcthrsh'] = float(config.get('avaProbs_thresholds_aps', 'lwcthrsh'))

    thresholds_aps['lwcthrsh_0'] = float(config.get('avaProbs_thresholds_aps', 'lwcthrsh_0'))
    thresholds_aps['lwcthrsh_1'] = float(config.get('avaProbs_thresholds_aps', 'lwcthrsh_1'))

    thresholds_aps['nsthrsh'] = float(config.get('avaProbs_thresholds_aps', 'nsthrsh'))
    thresholds_aps['release'] = release

    #thresholds for assigning avaprobs
    if scmopt == 'snp':
        threholds_ava = 'avaProbs_thresholds_avaprob_snp'
    elif scmopt == 'cro':
        threholds_ava = 'avaProbs_thresholds_avaprob_cro'
    else:
        raise ValueError('check opt in ini file, only snp or cro valid. your input is:', scmopt)
        
    thresholds_avaProb = {}
    thresholds_avaProb['damthrshnap'] = int(config.get(threholds_ava, 'damthrshnap'))
    thresholds_avaProb['inithrshnap'] = int(config.get(threholds_ava, 'inithrshnap'))
    thresholds_avaProb['propthrshnap'] = int(config.get(threholds_ava, 'propthrshnap'))
    thresholds_avaProb['damthrshpap'] = int(config.get(threholds_ava, 'damthrshpap'))
    
    thresholds_avaProb['inithrshpap'] = float(config.get(threholds_ava, 'inithrshpap'))
    thresholds_avaProb['propthrshpap'] = float(config.get(threholds_ava, 'propthrshpap'))
    thresholds_avaProb['drftthrsh'] = int(config.get(threholds_ava, 'drftthrsh'))
    thresholds_avaProb['dysisomax'] = int(config.get(threholds_ava, 'dysisomax'))
    thresholds_avaProb['lwcthrsh'] = float(config.get(threholds_ava, 'lwcthrsh'))

    thresholds_avaProb['lwcthrsh_0'] = float(config.get(threholds_ava, 'lwcthrsh_0'))
    thresholds_avaProb['lwcthrsh_1'] = float(config.get(threholds_ava, 'lwcthrsh_1'))
    thresholds_avaProb['release'] = release 

    #check if output dir is existing:
    isExist = os.path.exists(output_path)
    if not isExist:
        os.mkdir(output_path)
        print('output dir created...')
    else:
        print('output dir already existing...')

    # check if simulation folder already exists:
    isExist = os.path.exists(path)
    #create folder for output df:
    if not isExist:
        os.makedirs(path)
        print(' new directory - '+ station_name + ' -in output created...')
    else:
        print('dir for station-'+ station_name+ ' already existing...')

    #finding WL

    for ele in range(len(list_pro_red)):
        print('...algorithm startin for ', list_pro_red[ele],'...' )
        if rerun_find_WL == 1:
            print('...generate data from pro and smet file')
            # find weaklayers
            df_met, df_P = pro_find_aps.find_aps(season_start, season_end, drytime,wettime, thresholds_aps, list_pro_red[ele], list_smet_red[ele],read_in=rerun_find_WL, initilization_type=initilization_type)
            #save files as pkl
            

            df_P.to_pickle(path + '/df_P_'+ str(slope_aspects_sim[ele]) + '.pkl')
            df_met.to_pickle(path + '/df_met_' +str(slope_aspects_sim[ele]) + '.pkl')
        else:#load from pkl
            print('....load data from pkl files')
            df_P = pickle.load(open(path+'/df_P_'+ str(slope_aspects_sim[ele]) + '.pkl', "rb"))
            df_met = pickle.load(open(path +'/df_met_'+ str(slope_aspects_sim[ele]) + '.pkl', "rb"))


        #%% ------ from the output P (weak layer types, need instability THRSH to derive ava. problems!!) ------
    #% count cases: #days with dap (agingpersistentWL), #days with pap (persistent WL), #days with nap (non-persistent WL), date of wap (wet snow) onset
    #%tresholds from ROC curves!
    # assigne avaprobs to weaklayers

        if rerun_assign_avaprobs == 1:
            print('...generating avaprobs')
            df_P = pro_post_processing.assigne_avaprobs(df_P, thresholds_avaProb)
            df_P.to_pickle(path+ '/df_P_'+ str(slope_aspects_sim[ele])+'_avaprobs_'+'.pkl')
        else:
            print('... loading avaprobs')
            df_P = pickle.load(open(path +'/df_P_'+  str(slope_aspects_sim[ele]) + '_avaprobs_'+'.pkl', "rb"))
        
        print('... generartion of probs ' , str(slope_aspects_sim[ele]) ,'finished')
    print('...generations of avaprobs finished')
    # # graph de problems
    # start = datetime.datetime.strptime(season_start, '%Y-%m-%d')
    # end = datetime.datetime.strptime(season_end,'%Y-%m-%d')
    # color_napex = (120, 190,120)
    # color_napex = [element / 256 for element in color_napex]
    # color_winex = (0,100 ,20 )
    # color_winex = [element / 256 for element in color_winex]
    # color_papex = (65, 105 ,225 )
    # color_papex = [element / 256 for element in color_papex]
    # color_dapex = (25, 25, 112)
    # color_dapex = [element / 256 for element in color_dapex]
    # color_wapex = (150, 20, 0)
    # color_wapex = [element / 256 for element in color_wapex]
    # myLoc = mdates.MonthLocator()
    # myFmt = mdates.DateFormatter('%y-%b-%d') # %Y-%b-%d


    # #%% graph weather and problems and instability
    # fig, axs = plt.subplots(2,1, figsize = (15,6))
    # axs[0].semilogy(df_P['dy'], df_P['napDAM_extm2failMIN24_isrel'], 'x', color = color_napex, label= 'Dry non-persist.')
    # axs[0].semilogy(df_P['dy'], df_P['papDAM_extm2failMIN24_isrel'], 's', markerfacecolor='none', color= color_papex, label= 'Dry persist.')

    # axs[0].set_ylabel('Snow instability: t_f (h)')
    # axs[0].set_xlim([start.date(), end.date()])
    # axs[0].xaxis.set_major_formatter(myFmt)
    # axs[0].xaxis.set_major_locator(myLoc)
    # axs[0].set_xlim([start.date(), end.date()])
    # axs[0].legend(loc = 'lower left')

    # axs_0 = axs[0].twinx()

    # axs_0.plot(df_P['dy'], df_P['wapLWC_isrel']*100, 'o', color = color_wapex, markerfacecolor='none',label= 'Wet snow')
    # axs_0.legend(loc = 'lower left')

    # axs_0.set_ylabel('Snow instability: LWC_ind')
    # #axs[0].set_xlim([start.date(), end.date()])
    # #axs[0].xaxis.set_major_formatter(myFmt)
    # #axs[0].xaxis.set_major_locator(myLoc)
    # #axs[0].set_xlim([start.date(), end.date()])
    # axs_0.legend(loc = 'lower right')

    # #problems
    # axs[1].bar(df_P['dy'], df_P['napex_sele'], width = 0.4, bottom = 4, color = color_napex, label = 'New snow')
    # axs[1].bar(df_P['dy'], df_P['winex'], width = 0.4,bottom = 3, color = color_winex, label = 'Wind Slab')
    # axs[1].bar(df_P['dy'], df_P['papex_sele'],width = 0.4, bottom = 2, color = color_papex, label = 'Persistent' )
    # axs[1].bar(df_P['dy'], df_P['dapex_sele'], width = 0.4, bottom = 1, color = color_dapex, label ='Aging-persist.')
    # axs[1].bar(df_P['dy'], df_P['wapex_sele'], width = 0.4, bottom = 0, color = color_wapex, label ='Wet snow')

    # axs[1].set_ylabel('Ava problems')
    # axs[1].set_xlim([start.date(), end.date()])
    # axs[1].xaxis.set_major_formatter(myFmt)
    # axs[1].xaxis.set_major_locator(myLoc)
    # axs[1].set_xlim([start.date(), end.date()])
    # axs[1].legend(loc='lower left')

    # fig.tight_layout()
    # plt.savefig('avaprob.png')
    # plt.show()

    # # %% problem frequency

    # napno = np.sum(df_P['napex_sele'])
    # winno = np.sum(df_P['winex'])
    # papno = np.sum(df_P['papex_sele'])
    # dapno = np.sum(df_P['dapex_sele'])
    # wapno = np.sum(df_P['wapex_sele'])

    # PROBS = [napno, winno, papno, dapno, wapno]
    # labels = ['non-per.', 'wind-slab', 'per.' ,'aging-per', 'wet-snow']
    # colors = [color_napex, color_winex, color_papex, color_dapex, color_wapex]
    # fig1, ax1 = plt.subplots()
    # ax1.pie(PROBS/(sum(PROBS)/100),  labels = labels, colors =colors ,
    #         autopct='%1.1f%%', startangle=90, frame = True)
    # ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    # plt.savefig('avaprob_stat.png')
    # plt.show()

if __name__ == "__main__":
    """
    import debugpy
    debugpy.listen(5678)
    print('Waiting for debugger!')
    debugpy.wait_for_client()
    print('Attached!')
    """
    args = sys.argv[1:] # config_file = args[0]
    if (os.path.isfile(args[0])):
        main(args[0])
    else:
        name = str(args[0])
        print('File ({}) not found.'.format(name))