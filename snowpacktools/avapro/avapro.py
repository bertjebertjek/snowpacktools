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

# import matplotlib.pyplot as plt
# import matplotlib.dates as mdates
# from matplotlib.colors import BoundaryNorm, ListedColormap
# from matplotlib.ticker import AutoMinorLocator, FuncFormatter

from snowpacktools.avapro import find_aps
from snowpacktools.avapro import post_process_aps

def main(config_file):
    """Main function to run the algorithm on all .pro (+.smet) files
    in input folder defined in configuration file"""

    config = configparser.ConfigParser()
    config.read(config_file)

    """Parameters for research applications"""
    rerun_find_WL           = 1
    rerun_assign_avaprobs   = 1

    """Get list of available files"""
    SIM_FOLDER = config.get('AVAPRO', 'SIM_FOLDER_PATH')
    OUTPUT_FOLDER = os.path.join(SIM_FOLDER,'avapro-output/')
    list_pro  = sorted(glob.glob(SIM_FOLDER + "/*.pro"))
    list_smet = sorted(glob.glob(SIM_FOLDER + "/*.smet"))

    """Select certain aspects or filter for other things"""
    # SLOPE_ASPECTS=flat N E S W
    list_pro_red  = list_pro
    list_smet_red = list_smet

    if len(list_pro_red) == len(list_smet_red):
        pass
    elif len(list_pro_red) == 0:
        raise ValueError('[E]   No input, please check your input folder')
    else:
        raise ValueError('[E]   Number of .pro files and .smet files does not match')
        
     

    """Define output folder"""
    isExist = os.path.exists(OUTPUT_FOLDER)
    if not isExist:
        os.mkdir(OUTPUT_FOLDER)
        print('[I]  Output directory created')
    else:
        print('[I]  Output directory already exists')


    """Find potential weak layers (APS)"""
    print("[I]  Pro-list: ", list_pro_red)
    for ele in range(len(list_pro_red)):
        ele_name = list_pro_red[ele].split('/')[-1].split('.')[0]
        # print("[I]  Starting with PRO-file: ", list_pro_red[ele])
        print("[I]  Starting with PRO-file: ", ele_name)
        if rerun_find_WL == 1:
            print('[I]  Generate data from pro and smet files')
            df_met, df_P, meta_dict = find_aps.find_aps(config, list_pro_red[ele], list_smet_red[ele])
            
            """Save as pkl file"""
            df_P.to_pickle(os.path.join(OUTPUT_FOLDER,ele_name + '_df_P.pkl'))
            df_met.to_pickle(os.path.join(OUTPUT_FOLDER,ele_name + '_df_met.pkl'))
        else:
            print('[I]  Load data from pkl files')
            df_P = pickle.load(open(os.path.join(OUTPUT_FOLDER, ele_name + '_df_P.pkl'), "rb"))
            df_met = pickle.load(open(os.path.join(OUTPUT_FOLDER, ele_name + '_df_met.pkl'), "rb"))

        """Assign Avalanche Problems (APs)"""
        if rerun_assign_avaprobs == 1:
            print('[I]  Assigning avalanche problems from pkl files of tracked WLs')
            df_P = post_process_aps.assign_aps(df_P,config)
            df_P.to_pickle(os.path.join(OUTPUT_FOLDER, ele_name + '_df_P_APS.pkl'))
        else:
            print('[I]  Loading avalanche problems from pkl files')
            df_P = pickle.load(open(os.path.join(OUTPUT_FOLDER, ele_name + '_df_P_APS.pkl'), "rb"))
        
    print('[I]  Tracking WLs and assigning avalanche problems finished')


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