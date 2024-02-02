################################################################################
# Copyright 2024 Avalanche Warning Service Tyrol / Florian Herla               #
################################################################################
# This is free software you can redistribute/modify under the terms of the     #
# GNU Lesser General Public License 3 or later: http://www.gnu.org/licenses    #
################################################################################

import os
import sys
import configparser
import subprocess
import multiprocessing
import pkg_resources

import numpy as np
import pandas as pd

def aggregate(config):
    """Compute average/representative snow profiles from snowpack simulations stored in .pro files.

    The function can be used to aggregate profiles from an entire season, or it can be called day-by-day
    to aggregate the profiles as the eason proceeds in an operational mode.

    Parameters
    ----------
    config: configparser.ConfigParser
        A config as returned by the function `setup`. Can be modified by a custom config.ini file.
    
    Returns
    -------
    None:
        Instead of returning a python object, the function writes `.rds` and `.png` files stored in `./output/[...]`
    """

    """Create data frame that stores unique combinations of region, band, aspect
    and store into several csv files for parallel processing"""
    df = pd.read_csv(config.get('Paths', '_aggregates_vstations_csv_file'))
    dfuni = df[['region_id', 'band', 'aspect']].copy()
    dfuni = dfuni.drop_duplicates()
    dfuni_split = np.array_split(dfuni, config.getint('General','NTASKS'))
    for i in range(0,config.getint('General','NTASKS')):
        dfuni_split[i].to_csv(config.get('Paths','_aggregates_mp_csv') + str(i) + ".csv", index=False)

    """Run aggregation script in parallel"""
    print("[i]  Running aggregation script on multiple cpus.")
    print("")
    procs = []
    ## Start processes
    for i in range(0,config.getint('General','NTASKS')):
        proc = multiprocessing.Process(target=_worker_aggregation, args=(i, config))
        procs.append(proc)
        proc.start()
    ## Complete processes
    for proc in procs:
        proc.join()
    print("[i]  Number of cpus available: ", multiprocessing.cpu_count())
    print("[i]  Number of tasks used: ", config.getint('General','NTASKS'))
    print("[i]  Profile aggregation completed.")

    """Clean up"""
    for i in range(0,config.getint('General','NTASKS')):
        os.remove(config.get('Paths','_aggregates_mp_csv') + str(i) + ".csv")
    if config.getboolean('cleanup', 'dotinput'):
        os.rmdir("./input")



def setup(configfile, domain=''):
    """Set relevant config parameters to run aggregating function
    
    This function takes the default configuration of the package and updates the settings 
    with the provided configfile. It also creates the relevant directories if they don't exist yet.
    
    Parameters
    ----------
    configfile: str
        'path/to/config.ini'
    domain: str
        domain descriptor to be appended to directory names

    Returns
    -------
    config: an updated configparser.ConfigParser instance
    """
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(pkg_resources.resource_filename('snowpacktools', 'aggregatepro/aggregate.ini'))
    config.read(configfile)
    
    if len(domain) == 0:
        domain_appendix = ''
    else:
        domain_appendix = "-" + domain
    #domain_appendix = os.path.basename(config['Paths']['_snp_output_dir']).replace("snp", "")  # e.g. "-subtirol10000"
    
    ## The following ones could be made accessible to user config for more flexible control?!
    config['Paths']['_aggregates_output_dir'] = "./output/snp-aggregates" + domain_appendix
    config['Paths']['_aggregates_figures_dir'] = "./output/snp-aggregates-figs" + domain_appendix
    
    ## The following ones need no changing
    config['Paths']['_aggregates_mp_csv']         = './input/aggregates_mp'  # do not modify!
    config['Paths']['_aggregates_Rscript_path'] = pkg_resources.resource_filename('snowpacktools', 
                                                                                  'aggregatepro/aggregate_gridded_forecasts.R')
    config['Paths']['_aggregates_plotters_path'] = pkg_resources.resource_filename('snowpacktools', 
                                                                                   'aggregatepro/plotters.R')
    config['Paths']['_ini_runtime_domain'] = configfile  # already set in the context of 'awsome', but not for outside standalone use
    config['Aggregate']['DOMAIN'] = domain

    if config.get('Paths', '_aggregates_vstations_csv_file') == '_vstations_csv_file':
        config['Paths']['_aggregates_vstations_csv_file'] = config['Paths']['_vstations_csv_file']
    if config.get('Paths', '_aggregates_snp_pro_dir') == '_snp_ouput_dir':
        config['Paths']['_aggregates_snp_pro_dir'] = config['Paths']['_snp_ouput_dir']

    """Create directories"""
    os.makedirs("./output", exist_ok=True)
    os.makedirs(config['Paths']['_aggregates_output_dir'], exist_ok=True)
    os.makedirs(config['Paths']['_aggregates_figures_dir'], exist_ok=True)
    if os.path.exists("./input"):
        config.read_string("[cleanup]\ndotinput = False")
    else:
        os.makedirs("./input")
        config.read_string("[cleanup]\ndotinput = True")

    with open(configfile, "w") as cfgfile:
        config.write(cfgfile)

    return config



def _worker_aggregation(i, config):
    """Worker function that calls R script for aggregating gridded snow profiles stored in .pro files."""
    returnCode = subprocess.call(["Rscript", config.get('Paths', '_aggregates_Rscript_path'), 
                                  "--config", config.get("Paths","_ini_runtime_domain"), 
                                  "--mp_csv", config.get('Paths','_aggregates_mp_csv') + str(i) + ".csv"])
    if returnCode==0:
        print("[i]  Aggregation script successful for process number {}.".format(i))
    else:
        print("[E]  Aggregation script failed for process number {}!".format(i))



if __name__ == "__main__":
    """Manually aggregate gridded snow profiles from the command line.
    
    Usage
    -----
    python3 gridded.py configfile [domain]
        Make sure your working directory is set correctly before calling the module and your custom 
        configfile is adjusted to your environment.

    Returns
    -------
    The module will aggregate the profiles and store the results in subdirectories of './output/'.
    Directory names will contain the domain string.
    """

    if len(sys.argv) < 2 | len(sys.argv) > 3:
        sys.exit("[E] Synopsis: python3 gridded.py configfile [domain]")

    print(f'[i] Working directory set to {os.getcwd()}')
    
    configfile = sys.argv[1]
    # configfile = "/home/flo/documents/code/awsome/models/SNOWPACK/forecasts/input/forecast_runtime_subtirol10000.ini"

    if len(sys.argv) == 2:
        config = setup(configfile)
    elif len(sys.argv) == 3:
        domain = sys.argv[2]
        config = setup(configfile, domain)

    aggregate(config)