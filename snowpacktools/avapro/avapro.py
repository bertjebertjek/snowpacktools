################################################################################
# Copyright 2022 Avalanche Warning Service Tyrol                               #
################################################################################
# This is free software you can redistribute/modify under the terms of the     #
# GNU Lesser General Public License 3 or later: http://www.gnu.org/licenses    #
################################################################################

import os
import glob
import sys
import configparser
import pickle
from datetime import datetime, timedelta, date

from snowpacktools.avapro import find_aps, post_process_aps, visually_process_aps


def avapro(config_file):
    """Main function to run the algorithm on all .pro (+.smet) files
    in input folder defined in configuration file"""

    config = configparser.ConfigParser()
    config.read(config_file)

    """Get opera date"""
    datetime_format = '%Y-%m-%dT%Hh%M'
    # anatime = datetime.strptime('0600','%H%M').time()
    anatime = datetime.strptime('0000','%H%M').time()
    if config.get('AVAPRO','DATE_OPERA') == 'TODAY':
        datetime_today = datetime.combine(date.today(), anatime)
    else:
        date_today     = datetime.strptime(config.get('AVAPRO','DATE_OPERA'), "%Y-%m-%d").date()
        datetime_today = datetime.combine(date_today, anatime)

    """Make sure opera date is within season limits (otherwise print warning and set to end of season)"""
    date_season_end     = datetime.strptime(config.get('AVAPRO','SEASON_END'), "%Y-%m-%d").date()
    datetime_season_end = datetime.combine(date_season_end, anatime)
    if datetime_today > datetime_season_end:
        datetime_today = datetime_season_end
    datetime_today_str = datetime.strftime(datetime_today, datetime_format)

    """Parameters for research applications"""
    rerun_find_WL            = int(config["AVAPRO"]["rerun_find_WL"])
    rerun_assign_avaprobs    = int(config["AVAPRO"]["rerun_assign_avaprobs"])
    run_visualize_avaprobs   = int(config["AVAPRO"]["run_visualize_avaprobs"])

    """Output directories (one for pickle files, one for figures)"""
    OUTPUT_DIR      = config.get('AVAPRO', 'OUTPUT_DIR')
    OUTPUT_DIR_FIGS = config.get('AVAPRO', 'OUTPUT_DIR_FIGS')
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR_FIGS, exist_ok=True)

    """Input directories: Get list of available files"""
    if config.get('AVAPRO', 'SNP_FILE') == "NONE":
        list_pro  = sorted(glob.glob(config.get('AVAPRO', 'SNP_DIR') + "/*.pro"))
        list_smet = sorted(glob.glob(config.get('AVAPRO', 'SNP_DIR') + "/*.smet"))
        
        """Select certain aspects or filter for other things"""
        # SLOPE_ASPECTS=flat N E S W
    else:
        list_pro  = [os.path.join(config.get('AVAPRO', 'SNP_DIR'),config.get('AVAPRO', 'SNP_FILE'))]
        list_smet = [os.path.join(config.get('AVAPRO', 'SNP_DIR'),config.get('AVAPRO', 'SNP_FILE').split(".")[0] + ".smet")]
        
    if len(list_pro) == len(list_smet):
        pass
    elif len(list_pro) == 0:
        raise ValueError('[E]   No input, please check your input folder')
    else:
        raise ValueError('[E]   Number of .pro files and .smet files does not match')

    """Find potential weak layers (APS)"""
    print("[i]  Pro-list: ", list_pro)
    for ele in range(len(list_pro)):
        pro_name = list_pro[ele].split('/')[-1].split('.')[0]
        print("[i]  Running AVAPRO on PRO-file: ", pro_name)
        if rerun_find_WL == 1:
            print('[i]  Generate data from pro and smet files')
            df_met, df_P, meta_dict = find_aps.find_aps(config, list_pro[ele], list_smet[ele])
            
            """Save as pkl file"""
            df_P.to_pickle(os.path.join(OUTPUT_DIR,pro_name + '_df_P.pkl'))
            df_met.to_pickle(os.path.join(OUTPUT_DIR,pro_name + '_df_met.pkl'))
        else:
            print('[i]  Load data from pkl files')
            df_P = pickle.load(open(os.path.join(OUTPUT_DIR, pro_name + '_df_P.pkl'), "rb"))
            df_met = pickle.load(open(os.path.join(OUTPUT_DIR, pro_name + '_df_met.pkl'), "rb"))

        """Assign Avalanche Problems (APs)"""
        if rerun_assign_avaprobs == 1:
            print('[i]  Assigning avalanche problems from pkl files of tracked WLs')
            df_P = post_process_aps.assign_aps(df_P,config)
            df_P.to_pickle(os.path.join(OUTPUT_DIR, pro_name + '_df_P_APS.pkl'))
        else:
            print('[i]  Loading avalanche problems from pkl files')
            df_P = pickle.load(open(os.path.join(OUTPUT_DIR, pro_name + '_df_P_APS.pkl'), "rb"))
        
        """Visualize Avalanche Problems (APs)"""
        if run_visualize_avaprobs == 1:
            print('[i]  Visualizing avalanche problems from pkl files of tracked WLs')
            path_to_pro      = list_pro[ele]
            output_path      = os.path.join(OUTPUT_DIR_FIGS, pro_name + ".png")
            output_path_Punstable = os.path.join(OUTPUT_DIR_FIGS, pro_name + "_Punstable.png")
            output_path_sk38 = os.path.join(OUTPUT_DIR_FIGS, pro_name + "_sk38.png")
            res         = "1d"
            #visually_process_aps.plot_aps_and_profile_evolution(df_P, path_to_pro, DATETIME_STR=datetime_today_str, output_path=output_path, var='grain_type', res=res, second_var='Sk38', COLOR_SCHEME='SARPGR',DATE_RANGE=['NONE','NONE'])
            #visually_process_aps.plot_aps_and_profile_evolution(df_P, path_to_pro, DATETIME_STR=datetime_today_str, output_path=output_path, var='grain_type', res=res, second_var='NONE', COLOR_SCHEME='IACS2',DATE_RANGE=['NONE','NONE'])
            visually_process_aps.plot_aps_and_profile_evolution(df_P, path_to_pro, DATETIME_STR=datetime_today_str, output_path=output_path, var='grain_type', res=res, second_var='Punstable', COLOR_SCHEME='SARPGR',DATE_RANGE=['NONE','NONE'])
            # visually_process_aps.plot_aps_and_profile_evolution(df_P, path_to_pro, DATETIME_STR=datetime_today_str, output_path=output_path_sk38, var='grain_type', res=res, second_var='Sk38', COLOR_SCHEME='SARPGR',DATE_RANGE=['NONE','NONE'])

    print('[i]  Tracking WLs and assigning avalanche problems finished')


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
        avapro(args[0])
    else:
        name = str(args[0])
        print('File ({}) not found.'.format(name))