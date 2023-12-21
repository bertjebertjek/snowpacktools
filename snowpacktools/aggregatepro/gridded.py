#! /usr/bin/python3
################################################################################
# Copyright 2022 Avalanche Warning Service Tyrol                               #
################################################################################
# This is free software you can redistribute/modify under the terms of the     #
# GNU Lesser General Public License 3 or later: http://www.gnu.org/licenses    #
################################################################################

import os
import configparser
import subprocess
import multiprocessing
import pkg_resources

import numpy as np
import pandas as pd

def aggregate(config):

    """Update config"""
    ## The following ones should be made accessible to user config!
    domain_appendix = os.path.basename(config['Paths']['_snp_output_dir']).replace("snp", "")  # e.g. "-subtirol10000"
    config['Paths']['_aggregates_output_dir'] = "./output/snp-aggregates" + domain_appendix
    config['Paths']['_aggregates_figures_dir'] = "./output/snp-aggregates-figs" + domain_appendix
    
    ## The following ones need no changing
    config['Paths']['_aggregates_mp_csv']         = './input/aggregates_mp'  # do not modify!
    config['Paths']['_aggregates_Rscript_path'] = pkg_resources.resource_filename('snowpacktools', 
                                                                                  'aggregatepro/aggregate_gridded_forecasts.R')
    config['Paths']['_aggregates_plotters_path'] = pkg_resources.resource_filename('snowpacktools', 
                                                                                   'aggregatepro/plotters.R')

    if config['Paths']['_aggregates_ini'] == "NA":
         config['Paths']['_aggregates_ini'] = pkg_resources.resource_filename('snowpacktools', 'aggregatepro/aggregate.ini')
         
    with open(config.get("Paths","_ini_temporary"), 'w') as configfile:
        config.write(configfile)

    os.makedirs(config['Paths']['_aggregates_output_dir'], exist_ok=True)
    os.makedirs(config['Paths']['_aggregates_figures_dir'], exist_ok=True)
    if os.path.exists("./input"):
        del_dotinput = False
    else:
        os.makedirs("./input")

    """Create data frame that stores unique combinations of region, band, aspect
    and store into several csv files for parallel processing"""
    # miframe = pd.MultiIndex.from_product([
    #     df['region_id'].unique(), 
    #     df['band'].unique(), 
    #     df['aspect'].unique()
    # ], names=['region', 'band', 'aspect']).to_frame(index=False)
    # miframe_split = np.array_split(miframe, config.getint('Forecast','NTASKS'))
    # for i in range(0,config.getint('Forecast','NTASKS')):
    #     miframe_split[i].to_csv(config.get('Paths','_aggregates_mp_csv') + str(i) + ".csv", index=False)
    df = pd.read_csv(config.get('Aggregate', '_vstations_csv_file'))
    dfuni = df[['region_id', 'band', 'aspect']].copy()
    dfuni = dfuni.drop_duplicates()
    dfuni_split = np.array_split(dfuni, config.getint('Forecast','NTASKS'))
    for i in range(0,config.getint('Forecast','NTASKS')):
        dfuni_split[i].to_csv(config.get('Paths','_aggregates_mp_csv') + str(i) + ".csv", index=False)

    """Run aggregation script in parallel"""
    print("[i]  Running aggregation script on multiple cpus.")
    print("")
    procs = []
    ## Start processes
    for i in range(0,config.getint('Forecast','NTASKS')):
        proc = multiprocessing.Process(target=_worker_aggregation, args=(i, config))
        procs.append(proc)
        proc.start()
    ## Complete processes
    for proc in procs:
        proc.join()
    print("[i]  Number of cpus available: ", multiprocessing.cpu_count())
    print("[i]  Number of tasks used: ", config.getint('Forecast','NTASKS'))
    print("[i]  Profile aggregation completed.")

    for i in range(0,config.getint('Forecast','NTASKS')):
        os.remove(config.get('Paths','_aggregates_mp_csv') + str(i) + ".csv")
    if del_dotinput:
        os.rmdir("./input")



def _worker_aggregation(i, config):
    """Call R script for aggregating gridded snow profiles stored in .pro files."""
    rscript_path = pkg_resources.resource_filename('snowpacktools', 'aggregatepro/aggregate_gridded_forecasts.R')
    returnCode = subprocess.call(["Rscript", config.get('Paths', '_aggregates_Rscript_path'), 
                                  "--config", config.get("Paths","_ini_temporary"), 
                                  "--mp_csv", config.get('Paths','_aggregates_mp_csv') + str(i) + ".csv"])
    if returnCode==0:
        print("[i]  Aggregation script successful for process number {}.".format(i))
    else:
        print("[E]  Aggregation script failed for process number {}!".format(i))




if __name__ == "__main__":
    """
    Have to sort out config mess and declutter!
    """
    print("Not implemented yet")

    # config = configparser.ConfigParser()
    # config.read("/home/flo/documents/code/awsome/models/SNOWPACK/forecasts/input/forecast_localHF.ini")
    
    # aggregate(config)
    