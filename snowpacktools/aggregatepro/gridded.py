"""
This module is designed to aggregate simulated snow profiles for avalanche forecasting purposes.
It processes datasets of gridded snowpack simulations stored in .pro and .smet files and supports operations
both in a batch mode (or research mode) for past data and in an operational mode on a day-by-day 
basis during the snow season.

The module integrates functionalities from Python for data management and R for aggregating tasks
(making use of the R packages sarp.snowprofile and sarp.snowprofile.[-alignment & -pyface]),
leveraging multiprocessing to handle data across multiple CPU cores. The results of the
aggregation are saved in various formats including .rds for R data objects of representative profiles
and .png for visual plots.

Key functionalities include:
- Reading and processing custom configuration files for setting up the environment and settings.
- Filtering and grouping snow profile data sets based on geographical and temporal criteria.
- Running the data aggregation in parallel to optimize processing time.
- Writing the output directly to files.

Usage:
------
The module is optimized for usage with the [AWSOME](https://gitlab.com/avalanche-warning) framework, 
but can also be applied to other settings of gridded snowpack simulations. It can either be used as
a script with command line arguments, or can be run from within Python.

1. Script usage:
When calling the module as a script, command line arguments have to specify the config file and
(optional) domain name. Ensure the correct setup of the working environment (script assumes to be called 
from the correct directory, if file paths are given as relative and not absolute paths). Necessary 
file paths and settings can be set as outlined in the provided configuration file `aggregate.ini`.

Examples:
    $ python3 gridded.py config.ini
    $ python3 gridded.py config.ini domain_name (if domain name differs from config: [General] domain)

2. Interactive Python Usage:
Using the module interactively from within Python is very similar to calling it as a script, but it
does ship with a few more functionalities, such as aggregating custom polygons that differ from the 
regions specified in the vstations-csv file (see `aggregate.ini` for more information.)

Examples:
  1. General use
    >>> import configparser
    >>> from snowpacktools.aggregatepro import gridded as gag
    >>> config = configparser.ConfigParser()
    >>> config.read('your/custom/config.ini')
    >>> gag.aggregate(config)

  2. AWSOME use with existing domain
    >>> gag.aggregate(domain='domain')

  3. Custom polygons that overwrite regions in vstations-csv
    >>> gag.aggregate(config, group_regions_geojson='custom_polygons.geojson')

See documentation of `aggregate` for more info.

Dependencies:
-------------
  * Python libraries like pandas, geopandas, shapely, numpy, etc.
  * R packages 
    - (required): sarp.snowprofile, sarp.snowprofile.alignment, sarp.snowprofile.pyface, stringr, configr
    - (optional): progress
    - To make the `sarp.snowprofile.pyface` package work seemlessly, it is advised to define an enviroment 
      variable `RETICULATE_PYTHON` or `RETICULATE_PYTHON_ENV` that points to the python executable or the 
      python environment (venv or conda).
  * .pro and .smet files that store the snowpack simulations
    - make sure your version of SNOWPACK is not older than April 2024 to benefit from important changes in
      the AsciiIO write routine, see bullets below.
    - .smet file needs to contain ski_pen 
    - .pro file needs to contain deposition date information for each layer as code 0505. This can be achieved
      by the setting in the SNOWPACK config file: [Output] PROF_AGE_OR_DATE = DATE.
  * vstations-csv file that stores information about which profiles from which regions to be aggregated, e.g.
        vstation,easting,northing,lon,lat,elev,band,region_id
        100001A,615000,7728000,17.96,69.63,150,0000-0300,3011
        100002A,617000,7720000,18.00,69.56,75,0000-0300,3011
    `aspect` can also be included, but can also be handled automatically as long as the vstation is given as
    id ending on 'A' or '0' for flat fields. In that case, rows for the four cardinal aspects will be 
    generated automatically.
    
Output:
-------
The module will write files to `./input` and `./output` and also create the directories if they don't exist.

License:
--------
Distributed under the AGPLv3 license.

Authors:
--------
Florian Herla <fherla@sfu.ca>

Potential Future Improvements:
------------------------------
 * Currently, aggregation is only supported for daily sampling. Useful to adjust for any time sampling, or at least hourly, in the future.
 * When implementing hourly (or finer than daily) sampling, update 'valid' string in figures to datetime instead of date, and also include a 
   'computed' string (that represents DATE_OPERA).
 * Besides the static `.png` figures created by the current implementation, there is room for designing interactive visualizations that allow 
   forecasters to inspect various spatial distributions from the underlying individual profiles. Could be implemented.
     - To facilitate these interactions, layers from the average profile are backtracked to all underlying individual layers. 
       Therefore, all the data from the individual profiles needs to be stored in `.rds` files, which consumes a lot of space. 
       Either the underlying aggregation routine needs to be redesigned to allow for backtracking of layers to the original `.pro` files, 
       or the disk space needs to be available (at least for the current season).
"""

import os
import sys
import configparser
import subprocess
import multiprocessing
import pkg_resources
import shapely

import numpy as np
import pandas as pd
import geopandas as gpd

from datetime import datetime

import pdb


def aggregate(config=None, domain=None, group_regions_geojson=None, region_ids=None, aspects=None, bands=None):
    """
    Compute average/representative snow profiles from SNOWPACK simulations stored in .pro files.

    This function aggregates snowpack data from either a complete season or on a daily basis as needed 
    during the operational snow season. It's capable of customizing the aggregation based on geographic 
    and temporal parameters specified by the user. The function leverages multiprocessing to handle 
    large datasets efficiently across multiple CPUs. The results are stored in .rds and .png formats 
    within `./output/aggregates[-figs]-domain`.

    Parameters
    ----------
    config : configparser.ConfigParser, optional
        Configuration settings as returned by the function `setup`. If not provided, the configuration 
        is expected to be specified through the `domain` parameter.
    domain : str, optional
        A domain name for which the aggregation is to be performed. Within the context of AWSOME, the 
        function will find the relevant config files for existing domains if `config` is not provided. 
        The domain string is also used to name the output directories.
    group_regions_geojson : str, optional
        Path to a GeoJSON file defining custom polygons for which the aggregation should be performed.
        If provided, this overrides the regions specified in the vstations-csv file. Each polygon feature
        needs an `id` property. Ideally, the CRS is specified in the geojson, will be converted to WGS84.
    region_ids : list of int or str, optional
        Specific region IDs to limit the aggregation to certain geographic regions. Those are either the
        IDs from the vstation-csv file or from the GeoJSON.
    aspects : list of str, optional
        Specific aspects (e.g., flat, north, south, east, west) to filter the data aggregation.
    bands : list of str, optional
        Elevation bands to filter the data aggregation.

    Returns
    -------
    None
        The function does not return a value. Instead, it writes output directly to files at
        `./output/aggregates[-figs]-domain`, including both data files (.rds) and visual plots (.png).

    Raises
    ------
    FileNotFoundError
        If the configuration file cannot be found and no `config` is passed.

    Notes
    -----
    Ensure that the necessary files exist. The function also prints status messages to the console, 
    providing updates on the processing stages and any potential errors encountered during execution.
    For a documentation of the possible config settings, see `aggregate.ini`

    Example
    -------
    Running the function with custom configuration:

        >>> import configparser
        >>> config = configparser.ConfigParser()
        >>> config.read('path/to/custom/config.ini')
        >>> aggregate(config)

    Running with default settings for a given domain:

        >>> aggregate(domain='your_domain_name')

    Filtering specific locations:
        
        >>> aggregate(config, region_ids=['SnowyValley', 'AvyRidge'], bands=['treeline'], aspects=['north', 'south'])
        
        >>> aggregate(config, group_regions_geojson='custom_polygons.geojson', region_ids=['custom_region_id'])
    """


    if not config:
        configfile = os.path.expanduser('~snowpack/gridded-chain/input/runtime_') + domain + '.ini'
        if os.path.exists(configfile):
            config = setup(configfile)
        else:
            raise FileNotFoundError(f"The file '{configfile}' does not exist. Either make sure the file exists, or run the function with the config instead of the domain name, e.g. " +
                                    "\nconfig = gag.setup(configfile)  # your personal path/to/configfile.ini" +
                                    "\ngag.aggregate(config)")
    
    df_unique_grouping = determine_groupings(config, group_regions_geojson, region_ids, aspects, bands)
    df_unique_grouping.sort_values(['region_id', 'band', 'aspect']).to_csv(config.get('Paths','_aggregates_groupings_csv') + "summary.csv", index=False)
    print(df_unique_grouping)
    if df_unique_grouping.shape[0] == 0:
        raise ValueError("Cannot find a single grouping of vstations. Maybe all your vstation simulations were erroneous (i.e., vstations csv error column != 0)?")

    """set parallel processing context"""
    print(f"[i]  Running aggregation script on {config.getint('General','ncpus')} CPUs.\n\n")
    config['Aggregate']['ntasks'] = str(min(df_unique_grouping.shape[0], config.getint('General', 'ntasks')))
    config['Aggregate']['ncpus_per_task'] = str(max(1, config.getint('General', 'ncpus') // config.getint('Aggregate', 'ntasks')))
    with open(config.get('Paths', '_ini_runtime_domain'), "w") as cfgfile:
        config.write(cfgfile)
    ## Store one csv file with groupings for each task to process
    df_unique_grouping_split = np.array_split(df_unique_grouping, config.getint('Aggregate','ntasks'))
    for i in range(0,config.getint('Aggregate','ntasks')):
        df_unique_grouping_split[i].to_csv(config.get('Paths','_aggregates_groupings_csv') + str(i) + ".csv", index=False)

    """Run aggregation script in parallel"""
    procs = []
    t0 = datetime.now()
    ## Start processes
    for i in range(0,config.getint('Aggregate','ntasks')):
        # proc = _worker_aggregation(i, config)  # for easier debugging..
        proc = multiprocessing.Process(target=_worker_aggregation, args=(i, config))
        procs.append(proc)
        proc.start()
    ## Complete processes
    for proc in procs:
        proc.join()
    duration = datetime.now() - t0
    print("\n[i]  Number of cpus available: ", multiprocessing.cpu_count())
    print("[i]  Number of cpus used: ", config.getint('General','ncpus'))
    print("[i]  Number of tasks distributed: ", config.getint('Aggregate','ntasks'))
    # print("[i]  Number of cpus-per-task: ", config.getint('Aggregate','ncpus_per_task'))  # not implemented yet in Rscript, but can be!
    print(f"[i]  Profile aggregation completed in {':'.join(str(duration).split(':')[:2])} hours.")

    """Clean up"""
    for i in range(0,config.getint('Aggregate','ntasks')):
        os.remove(config.get('Paths','_aggregates_groupings_csv') + str(i) + ".csv")
    if config.getboolean('cleanup', 'dotinput'):
        os.rmdir("./input")
    if '_aggregates_vstations_csv_file_processed' in config.options('cleanup'):
        if config.getboolean('cleanup', '_aggregates_vstations_csv_file_processed'):
            os.remove(config.get('Paths', '_aggregates_vstations_csv_file_processed'))
            if config.get('Paths', '_aggregates_vstations_csv_file_processed').endswith('.full'):
                config['Paths']['_aggregates_vstations_csv_file_processed'] = config.get('Paths', '_aggregates_vstations_csv_file_processed')[:-5]
                with open(config.get('Paths', '_ini_runtime_domain'), "w") as cfgfile:
                    config.write(cfgfile)


def determine_groupings(config, group_regions_geojson, region_ids, aspects, bands):
    """
    Filters and prepares grouping information for snow profile aggregation.

    Parameters
    ----------
    See `aggregate` for a detailed description.

    Returns
    -------
    pandas.DataFrame
        DataFrame with unique combinations of region, band, and aspect to be aggregated.
    """

    df = pd.read_csv(config.get('Paths', '_aggregates_vstations_csv_file_processed'))

    """Generate new vstations_csv in case a group_regions_geojson is provided"""
    if group_regions_geojson:
        df['geometry'] = df.apply(lambda row: shapely.geometry.Point(row['lon'], row['lat']), axis=1)
        gdf = gpd.GeoDataFrame(df, geometry='geometry', crs="EPSG:4326") # set epsg code for WGS84 format
        polygons = gpd.read_file(group_regions_geojson)
        polygons = polygons.to_crs(epsg=4326)
        gdf = gpd.sjoin(gdf, polygons, how='inner', predicate='within')
        gdf['region_id'] = gdf['id']
        df = gdf[[col for col in ['vstation', 'easting', 'northing', 'lon', 'lat', 'elev', 'band', 'region_id', 'error'] if col in gdf.columns]]
        # handle file extensions (`-ext` for external groupings, `.csv` default file extension, `.full` for aspect-laden files, e.g.:
        #                         /path/to/basename-ext.csv.full)
        root_vstationsfile, ext_vstationsfile = os.path.splitext(config.get('Paths', '_aggregates_vstations_csv_file_processed'))
        if ext_vstationsfile == '.full':
            root_vstationsfile, ext_prev_vstationsfile = os.path.splitext(root_vstationsfile)
            ext_vstationsfile = ext_prev_vstationsfile + ext_vstationsfile
        if not root_vstationsfile[-4:] == '-ext':
            config['Paths']['_aggregates_vstations_csv_file_processed'] = root_vstationsfile + '-ext' + ext_vstationsfile

    """Extend DataFrame with aspect info if not present"""
    if 'aspect' not in df.columns:
        if config.get('General', 'nslopes') == '5':
            aspect_map = {
            'A': 'flat',
            '0': 'flat',
            '1': 'north',
            '2': 'east',
            '3': 'south',
            '4': 'west'
            }
            aspect_suffix = pd.Series(['', 1, 2, 3, 4]*len(df))
            df = pd.concat([df]*5).sort_index(kind='merge').reset_index(drop=True)
            df['vstation'] += aspect_suffix.astype(str)
            df['aspect'] = df['vstation'].str[-1].map(aspect_map)

            config['Paths']['_aggregates_vstations_csv_file_processed'] = config.get('Paths', '_aggregates_vstations_csv_file_processed') + '.full'
            config['cleanup']['_aggregates_vstations_csv_file_processed'] = 'True'
        else:
            raise NotImplementedError("Nslopes other than 5 are not implemented yet.")
        
    if 'error' in df.columns:
        df_err = df.loc[df['error'] > 0, ].copy()
        df = df.loc[df['error'].isin([0]), ].copy()

    with open(config.get('Paths', '_ini_runtime_domain'), "w") as cfgfile:
            config.write(cfgfile)
    df.to_csv(config.get('Paths', '_aggregates_vstations_csv_file_processed'), index=False)

    """Filter region_ids, bands, aspects"""
    if region_ids:
        df = df[df['region_id'].isin(region_ids)]
    if aspects:
        df = df[df['aspect'].isin(aspects)]
    if bands:
        df = df[df['band'].isin(bands)]

    """Create groupings csv that stores unique combinations of region, band, aspect"""
    dfuni = df[['region_id', 'band', 'aspect']].copy()
    dfuni = dfuni.drop_duplicates()

    ## arrange the data frame in a way that alternates combinations with few and many stations, 
    ## so that resulting workers have roughly the same amount of work
    station_counts = df.groupby(['region_id', 'band', 'aspect'])['vstation'].nunique().reset_index(name='nstations')
    if 'df_err' in locals():
        err_stations= df_err.groupby(['region_id', 'band', 'aspect'])['error'].sum().reset_index(name='excl_nstations_err')
        station_counts = pd.merge(station_counts, err_stations, on=['region_id', 'band', 'aspect'], how='left')
        station_counts['excl_nstations_err'] = station_counts['excl_nstations_err'].fillna(0).astype(int)
    dfuni = pd.merge(dfuni, station_counts, on=['region_id', 'band', 'aspect'], how='left')
    dfuni = dfuni.sort_values(by='nstations').reset_index(drop=True)
    dfuni = dfuni.astype(str)
    midpoint = len(dfuni) // 2
    if len(dfuni) % 2 == 0:
        # Even number of rows: interleave directly
        dfuni_alternating = pd.concat([dfuni[:midpoint].reset_index(drop=True),
                                    dfuni[midpoint:].reset_index(drop=True).iloc[::-1]], axis=1).stack().reset_index(drop=True)
    else:
        # Odd number of rows: handle the last row separately
        dfuni_alternating = pd.concat([dfuni[:midpoint+1].reset_index(drop=True),
                                    dfuni[midpoint+1:].reset_index(drop=True).iloc[::-1]], axis=1).stack().reset_index(drop=True)
        
    dfuni = pd.DataFrame(dfuni_alternating.values.reshape(dfuni.shape), columns = dfuni.columns)
    
    return dfuni



def setup(configfile, domain=''):
    """
    Initializes and updates configuration settings from a file for the aggregation process.

    The user provided config file can only contain custom modifications. All missing settings
    will be taken from the template file containing all defaults, `aggregate.ini`.

    Parameters
    ----------
    configfile : str
        File path to the configuration file.
    domain : str, optional
        Domain descriptor to be appended to directory names (defaults to value in config: [General] domain).

    Returns
    -------
    configparser.ConfigParser
        An updated configuration object reflecting the settings from the provided file.
    """
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(pkg_resources.resource_filename('snowpacktools', 'aggregatepro/aggregate.ini'))
    config.read(configfile)
    
    if len(domain) == 0:
        domain = config.get('General', 'domain')
    domain_appendix = "-" + domain
    
    ## The following settings need no changing
    config['Paths']['_aggregates_groupings_csv']         = f'./input/aggregates_groupings/{domain}-'  # do not modify!
    config['Paths']['_aggregates_rscript_path'] = pkg_resources.resource_filename('snowpacktools', 
                                                                                  'aggregatepro/aggregate_gridded_forecasts.R')
    config['Paths']['_aggregates_plotters_path'] = pkg_resources.resource_filename('snowpacktools', 
                                                                                   'aggregatepro/plotters.R')
    config['Paths']['_ini_runtime_domain'] = configfile  # already set in the context of 'awsome', but not for outside standalone use

    if config.get('Paths', '_aggregates_output_basedir') == '_output_dir':
        try:
            config['Paths']['_aggregates_output_basedir'] = config.get('Paths', '_output_dir')
        except KeyError:
            config['Paths']['_aggregates_output_basedir'] = "./output"

    config['Paths']['_aggregates_output_dir'] = config['Paths']['_aggregates_output_basedir'] + "/aggregates" + domain_appendix
    config['Paths']['_aggregates_figures_dir'] = config['Paths']['_aggregates_output_basedir'] + "/aggregates-figs" + domain_appendix
    
    if config.get('Paths', '_aggregates_vstations_csv_file') == '_vstations_csv_file':
        config['Paths']['_aggregates_vstations_csv_file'] = config['Paths']['_vstations_csv_file']
    config['Paths']['_aggregates_vstations_csv_file_processed'] = config['Paths']['_aggregates_vstations_csv_file']
    if config.get('Paths', '_aggregates_snp_pro_dir') == '_snp_output_dir':
        config['Paths']['_aggregates_snp_pro_dir'] = config['Paths']['_snp_output_dir']

    """Create directories"""
    # os.makedirs("./output", exist_ok=True)
    os.makedirs(config['Paths']['_aggregates_output_dir'], exist_ok=True)
    os.makedirs(config['Paths']['_aggregates_figures_dir'], exist_ok=True)
    if os.path.exists("./input"):
        config.read_string("[cleanup]\ndotinput = False")
    else:
        os.makedirs("./input")
        config.read_string("[cleanup]\ndotinput = True")
    os.makedirs("./input/aggregates_groupings", exist_ok=True)

    with open(configfile, "w") as cfgfile:
        config.write(cfgfile)

    return config


def _worker_aggregation(i, config):
    """Worker function that calls R script for aggregating gridded snow profiles stored in .pro files."""
    returnCode = subprocess.call(["Rscript", config.get('Paths', '_aggregates_rscript_path'), 
                                "--config", config.get("Paths","_ini_runtime_domain"), 
                                "--mp_csv", config.get('Paths','_aggregates_groupings_csv') + str(i) + ".csv",
                                "--worker_int", str(i)])
    if returnCode==0:
        print("[i]  Aggregation script successful for process number {}.".format(i))
    else:
        print("[E]  Aggregation script failed at least partially for process number {}!".format(i))



if __name__ == "__main__":

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