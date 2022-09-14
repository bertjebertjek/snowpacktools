################################################################################
# Copyright 2022 Avalanche Warning Service Tyrol                               #
################################################################################
# This is free software you can redistribute/modify under the terms of the     #
# GNU Lesser General Public License 3 or later: http://www.gnu.org/licenses    #
################################################################################

import os
import sys
import time
import configparser

import numpy as np
import pandas as pd
import xarray # needed for time axis

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.ticker import AutoMinorLocator, FuncFormatter

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

import awset
import awsmet
from snowpacktools.snowpro import pro_helper

def read_pro(path,res='1h'):
    """Reads a .PRO file and returns a list of dataframes representing the evolving state of the snowpack.
    
    Arguments:
        path (str):             String pointing to the location of the .PRO file to be read
    Returns:
        snowpro_list (list):    List of dfs with each df representing one snow profile (one timestamp), column=layer, row=variables
        meta_dict:              Dictionary with metadata of snow profile
    """
    w, hours = pro_helper.set_resolution(res)

    PRO_CODE_DICT, VAR_CODES = pro_helper.get_pro_code_dict()
    variables = {}
    for var in VAR_CODES:
        variables[PRO_CODE_DICT[var]] = []
    
    meta_dict = {}
    
    """Open the PRO file and generate dict of variables with list of lines for each variable"""
    start_read_file = time.time()
    with open(path, "r") as f:
        file_content = f.readlines()
    
    section = '[STATION_PARAMETERS]'
    for line in file_content:
        line = line.rstrip('\n')

        if line == '[HEADER]':
            section = '[HEADER]'
            continue
        elif line == '[DATA]':
            section = '[DATA]'
            timestamp_of_interest = False

            """Drop variables that are not present in header of .PRO file"""
            keys_to_drop = []
            for key in variables:
                if len(variables[key]) == 0:
                    print('Variable ', key, ' is not found in the header of this .PRO file. It is dropped.')
                    keys_to_drop.append(key)
            
            for key in keys_to_drop:
                variables.pop(key)
            continue
        
        if section=='[STATION_PARAMETERS]':
            line_arr = line.split('= ')
            if len(line_arr) == 2:
                meta_dict[line_arr[0]] = line_arr[1]
        
        elif section=='[HEADER]':
            if line[:4] in VAR_CODES:
                variables[PRO_CODE_DICT[line[:4]]].append(line)

        else:
            if line[:4] == '0500':
                """Check that all variable lists are same length..."""
                n = len(variables['date'])
                for key in variables:
                    if (n-len(variables[key]))==1:
                        variables[key].append('-999')

                """Check if timestamp is of interest"""
                if int(line[-8:-6]) in hours:
                    timestamp_of_interest=True
                    # Add line to variable
                    variables['date'].append(line)
                else:
                    timestamp_of_interest=False

            elif line[:4] in VAR_CODES and timestamp_of_interest:
                variables[PRO_CODE_DICT[line[:4]]].append(line)

    """Check again that all variable lists are same length... (for no snow at end of season)"""
    n = len(variables['date'])
    for key in variables:
        if (n-len(variables[key]))==1:
            variables[key].append('-999')

    end_read_file = time.time()
    print('Reading of lines took: {}s'.format(int(end_read_file-start_read_file)))

    """Remove the header data (leave this, because it covers wrong user input with var_codes)"""
    for variable in variables.keys():
        try:
            variables[variable].pop(0)
        except:
            print('Attention: No values for', variable,'in your .pro file. Remove it out of var_code dictionary')

    """Check existence of soil layers (!exist for negative height values!)"""
    line_series = variables['height_m'][0].split(",")
    nvars = int(line_series[1])
    soil_vars     = []
    i_ground_surf = 0
    if nvars > 1:
        datapoints = line_series[-nvars:]
        datapoints = list(map(float, datapoints))
        if datapoints[0] < 0:
            print('Soil layers detected. The following variables contain values for soil:')
            i_ground_surf = datapoints.index(0.0)
            for varname in variables.keys():
                if varname != 'date':
                    line_series = variables[varname][0].split(",")
                    n = int(line_series[1])
                    if n+1==nvars:
                        soil_vars.append(varname)
            print(soil_vars)
    
    """Generate snow profile dataframe for each timestamp"""
    start_processing = time.time()
    snowpro_list = [pro_helper.snowpro_from_snapshot(i, variables, i_ground_surf, soil_vars) for i in range(len(variables['date']))]

    """Transform SLF graintype code into ICSSG standard abbreviation"""
    for df in snowpro_list:
        df['graintype'] = df['grain type (Swiss Code F1F2F3)'].apply(pro_helper.slf_graintype_to_ICSSG)

        i = np.arange(0,len(df)-1)
        df['bottom'] = 0
        df.loc[i,'bottom'] = df.loc[i+1,'height_m'].values
    end_processing = time.time()
    print('Generation of dataframes took: {}s'.format(int(end_processing-start_processing)))

    return snowpro_list, meta_dict


def get_smet_df(path):
    """Generates dataframe for further processing out of .smet file."""

    var = pro_helper.get_var_smet(path)

    df_smet = pd.read_csv(path, sep=" ", skiprows=18, skipinitialspace=True, names =var) 
    df_smet = df_smet.replace(-999.0, np.NaN)

    """Reduce dataframe to variables of interest"""
    variables_of_intrest = ['timestamp', 'TSS_mod', 'TSS_meas', 'T_bottom' ,'TSG','VW','DW','wind_trans24','VW_drift', 'MS_Wind', 'HS_mod', 'HS_meas',
                        'MS_Snow','hoar_size', 'HN72_24','HN24','MS_Rain','SWE','MS_Water']
    df_smet = df_smet.loc[:, df_smet.columns.intersection(variables_of_intrest)]

    df_smet['HS_mod']  = df_smet['HS_mod']/100
    df_smet['HS_meas'] = df_smet['HS_meas']/100
    df_smet['HN72_24'] = df_smet['HN72_24']/100
    df_smet['HN24']    = df_smet['HN24']/100
    return df_smet
    

def plot_single_profile(path_to_pro, DATETIME_STR,output_dir='output/', COLOR_SCHEME='IACS2',ax=None):
    """Plots snow profile of PRO file @DATETIME (shows hardness profile + grain type).
    
    Arguments:
        path_to_pro: path to PRO file to be plotted
    Returns:
        Figure
    """

    start_readin = time.time()
    df_pro_list_temp, meta_dict = read_pro(path_to_pro)
    end_readin = time.time()
    print('Reading of PRO file completed in {}s'.format(int(end_readin-start_readin)))

    LABELS_GRAIN_TYPE, COLORS_GRAIN_TYPE, HATCHES_GRAIN_TYPE, LABELS_GRAIN_TYPE_BAR, COLORS_GRAIN_TYPE_BAR, HATCHES_GRAIN_TYPE_BAR = pro_helper.get_grain_type_colors(COLOR_SCHEME)

    """Get closest profile of PRO file to DATETIME"""
    DATETIME = pd.Timestamp(DATETIME_STR)
    for df in df_pro_list_temp:
        if df.date.iloc[0] == DATETIME:
            print('Snow profile found for', DATETIME)
            # - Use df from now on - # 
            break

    hand_hardness_dict = pro_helper.get_hand_hardness_N_dict()
    df['hand_hardness_N'] = df['hand hardness']
    ZERO_HH_VAL = 50

    if ax==None:
        fig, ax = plt.subplots(1,1,figsize=(6,7))
    else:
        fig=None

    col_dict_labels     = dict(zip(LABELS_GRAIN_TYPE, COLORS_GRAIN_TYPE))
    hatches_dict_labels = dict(zip(LABELS_GRAIN_TYPE, HATCHES_GRAIN_TYPE))

    n_bar = len(LABELS_GRAIN_TYPE_BAR)
    col_nums = np.arange(0,n_bar)
    col_dict = dict(zip(col_nums, COLORS_GRAIN_TYPE_BAR[::-1]))
    cmap = ListedColormap([col_dict[x] for x in col_dict.keys()])

    cols=[]
    hatches=[]
    for row in range(0,len(df)):
        cols.append(col_dict_labels[df.loc[row,'graintype'][0]])
        hatches.append(hatches_dict_labels[df.loc[row,'graintype'][0]])
        df.loc[row,'hand_hardness_N'] = hand_hardness_dict[df.loc[row,'hand hardness']]
    bar_plot   = ax.barh(df['bottom'], df['hand_hardness_N'], height=df['thickness_m'], align='edge', color=cols, hatch=hatches) # label=labels[i])
    bar_plot_p = ax.barh(df['bottom'], ZERO_HH_VAL, height=df['thickness_m'], align='edge', color=cols, hatch=hatches)
        

    """COLORBAR (Norm, bins, formatter, ticks - lots of stuff to make colorbar look nice)"""
    if fig!=None:
        lulu = np.zeros((n_bar,n_bar))
        for nn,k in enumerate(col_dict.keys()):
            lulu[nn, :] = np.nan
        norm_bins = np.sort([*col_dict.keys()]) + 0.5
        norm_bins = np.insert(norm_bins, 0, np.min(norm_bins) - 1.0)

        norm = BoundaryNorm(norm_bins, n_bar, clip=True)
        fmt = FuncFormatter(lambda x, pos: LABELS_GRAIN_TYPE_BAR[::-1][norm(x)])
        diff = norm_bins[1:] - norm_bins[:-1]
        tickz = norm_bins[:-1] + diff / 2

        # contf = ax.contourf(lulu,cmap=cmap,norm=norm,levels=norm_bins) # just for colorbar
        contf = ax.contourf(lulu,cmap=cmap,norm=norm,levels=norm_bins,hatches=HATCHES_GRAIN_TYPE_BAR[::-1]) # just for colorbar
        cbar = fig.colorbar(contf, format=fmt, ticks=tickz,location='left', pad=0.04) # shrink=0.7, ax=[axes[1],axes[3], axes[5]]
        cbar.ax.grid(visible=False)

    """Temperature axis"""
    ax_t = ax.twiny()
    ax_t.plot(df['temperature'], df['height_m'], color='#DC143C',lw=1.5)
    ax_t.grid(visible=False)
    ax_t.xaxis.tick_bottom()
    ax_t.xaxis.set_label_position("bottom")
    ax_t.set_xlabel("snow temperature / °C", color='#DC143C')
    ax_t.tick_params(axis='x', colors='#DC143C')
    ax_t.set_xlim(min(df['temperature'])-2,0)

    """AXES"""
    XLIM = [-1100,50] 
    ax.set_xlim(XLIM)
    if fig!=None:
        ax.set_ylim(0,df.loc[0,'height_m']+0.1)
        ax.set_ylabel("height / m")
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
        
    # fig.canvas.draw()
    # tickz = ax.get_xticks()
    # tick_labels = ax.get_xticklabels()
    tickz       = [-1000,-750,-500,-250,0]
    tick_labels = [1000,750,500,250,0]
    ax.set_xticks(tickz)
    ax.set_xticklabels(tick_labels)
    
    ax_hh = ax.twiny()
    ax_hh.set_xlim(XLIM)
    ax_hh.grid(visible=False)
    ax_hh.xaxis.tick_top()
    ax_hh.xaxis.set_label_position("top")
    tickz       = [-1000,-500,-250,-100,-20]
    tick_labels = ['K','P','1F','4F','F']
    ax_hh.set_xticks(tickz)
    ax_hh.set_xticklabels(tick_labels)
    ax_hh.tick_params(axis="x",direction="in", pad=-18)

    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.set_xlabel("hand hardness / N")

    ax.yaxis.set_minor_locator(AutoMinorLocator())

    """Include Meta data in top left corner and save figure"""
    if fig!=None:
        header_str = 'Location:      ' + meta_dict['StationName'] + ' (' + DATETIME_STR + ')\nElevation:     ' + meta_dict['Altitude'] + \
                    'm\nSlope Angle: ' + str(int(float(meta_dict['SlopeAngle']))) + '°\nAspect:         ' + str(int(float(meta_dict['SlopeAzi'])))  + '°'
        ax.text(0.04,0.95,header_str,horizontalalignment='left',
                verticalalignment='top', fontsize=10,transform=ax.transAxes) # ma='left'

        # --- Save figure --- #
        fig_title = 'snow-profile-' + meta_dict['StationName'] + '-' +  DATETIME_STR + '.png'
        fig.tight_layout()
        fig.savefig(os.path.join(output_dir,fig_title), facecolor='w', edgecolor='w',
                    format='png', dpi=150)

        end_plotting = time.time()
        print('Visualization of snowpack evolution completed in {}s'.format(int(end_plotting-end_readin)))
    else:
        ax.text(0.04,0.93,DATETIME_STR,horizontalalignment='left',
                verticalalignment='top', fontsize=10, transform=ax.transAxes) # ma='left'
        return ax


def plot_snp_evo(path_to_pro, output_dir='output/', DATETIME_STR=None, var='grain_type', res='1h', second_var='NONE', COLOR_SCHEME='IACS2', DATE_RANGE=['NONE','NONE']):
    """Plots snowpack evolution (PRO-file). Different variables or grain type can be visualized and overlayed.
    
    Arguments:
        path_to_pro (str):    Path to PRO file to be plotted
        output_dir (str):     
        DATETIME_STR (str):
        var: optional   Can be used for visualizing other variables than grain type
        res (str):
        second_var (str):
        COLOR_SCHEME (str):
        DATE_RANGE (list):
    Returns:
        Figure
    """

    start_readin = time.time()
    df_pro_list_temp, meta_dict = read_pro(path_to_pro,res=res)
    end_readin = time.time()
    print('Reading of PRO file completed in {}s'.format(int(end_readin-start_readin)))

    """Filter for certain resolution and time frame"""
    w, hours = pro_helper.set_resolution(res)
    LABELS_GRAIN_TYPE, COLORS_GRAIN_TYPE, HATCHES_GRAIN_TYPE, LABELS_GRAIN_TYPE_BAR, COLORS_GRAIN_TYPE_BAR, HATCHES_GRAIN_TYPE_BAR = pro_helper.get_grain_type_colors(COLOR_SCHEME)
    RANGE_DICT = pro_helper.get_range_dict()

    df_pro_list = []
    for df in df_pro_list_temp:
        if df.date.iloc[0].hour in hours:
        # if (df.iloc[0].date.strftime('%m.%d') > season_start) or (df.iloc[0].dates.strftime('%m.%d') < season_end):
            df_pro_list.append(df)

    """COLOR MAP AND PREPROCESSING"""
    if var=='grain_type':
        col_dict_labels     = dict(zip(LABELS_GRAIN_TYPE, COLORS_GRAIN_TYPE))
        hatches_dict_labels = dict(zip(LABELS_GRAIN_TYPE, HATCHES_GRAIN_TYPE))

        n_bar = len(LABELS_GRAIN_TYPE_BAR)
        col_nums = np.arange(0,n_bar)
        col_dict = dict(zip(col_nums, COLORS_GRAIN_TYPE_BAR[::-1]))
        cmap = ListedColormap([col_dict[x] for x in col_dict.keys()])
    
    else:
        if var in ['Sk38','Sn38']:
            cmap_var = plt.get_cmap('plasma')
            # cmap_var = plt.get_cmap('BuPu_r')
        else:
            # cmap_var = plt.get_cmap('BuPu')
            cmap_var = plt.get_cmap('plasma_r')
        clev_var = np.linspace(RANGE_DICT[var][0],RANGE_DICT[var][1],11)
        cnorm_var = BoundaryNorm(boundaries=clev_var, ncolors=cmap_var.N, clip=True)

    var_alpha=1
    if second_var!='NONE':#
        if second_var in ['RTA']:
            cmap_var2 = pro_helper.get_whiteout_cmap()
            var_alpha=0.5
        else:
            # cmap_var2  = plt.get_cmap('gist_gray')
            cmap_var2 = pro_helper.get_whiteout_cmap(reverse=True)
            var_alpha=0.5
        clev_var2 = np.linspace(RANGE_DICT[second_var][0],RANGE_DICT[second_var][1],11)
        cnorm_var2 = BoundaryNorm(boundaries=clev_var2, ncolors=cmap_var2.N, clip=False)

    """VISUALIZATION"""
    if DATETIME_STR==None:
        fig, ax = plt.subplots(1,1,figsize=(14,6))
    else:
        fig, (ax,ax_prof) = plt.subplots(1,2,figsize=(14,6),sharey=True,gridspec_kw={'width_ratios': [3, 1]})

    h_max = []
    dates = []
    for df in df_pro_list:
        h_max.append(df.loc[0,'height_m'])
        dates.append(df.loc[0,'date'])
            
        if var=='grain_type':
            cols=[]
            hatches=[]
            for row in range(0,len(df)):
                cols.append(col_dict_labels[df.loc[row,'graintype'][0]])
                hatches.append(hatches_dict_labels[df.loc[row,'graintype'][0]])
            bar_plot = ax.bar(df.date[0], df['thickness_m'], width=w, bottom=df['bottom'], align='edge', color=cols, hatch=hatches, alpha=1) # label=labels[i])
        else:
            var_data = (df[var].values-RANGE_DICT[var][0])/(RANGE_DICT[var][1]-RANGE_DICT[var][0])
            cols = cmap_var(var_data)
            bar_plot = ax.bar(df.date[0], df['thickness_m'], width=w, bottom=df['bottom'], align='edge', color=cols)
        
        if second_var!='NONE':
            var_data2 = (df[second_var].values-RANGE_DICT[second_var][0])/(RANGE_DICT[second_var][1]-RANGE_DICT[second_var][0])
            cols2 = cmap_var2(var_data2)
            ax.bar(df.date[0], df['thickness_m'], width=w, bottom=df['bottom'], align='edge', color=cols2, alpha=0.8)

    """Line along snow surface"""
    if second_var!='NONE':
        # h_max = np.where(h_max == np.nan, 0, h_max)
        ax.plot(dates,h_max,ds='steps-post',lw=0.8,color='black', ls='--', alpha=0.67)

    """Hardness profile to the right"""
    if DATETIME_STR!=None:
        ax_prof = plot_single_profile(path_to_pro,DATETIME_STR,COLOR_SCHEME=COLOR_SCHEME,ax=ax_prof)
        ax.axvline(x=DATETIME_STR,ymin=-0.1, ymax=1.1, color='black', lw=3, ls='--')

    """COLORBAR (Norm, bins, formatter, ticks - lots of stuff to make colorbar look nice)"""
    if second_var!='NONE':
         # - Colorbar for second layer - # 
        n_var = 9
        lulu = np.zeros((n_var,n_var))
        for nn in range(0,n_var):
            lulu[nn, :] = np.nan # nn
        contf = ax.contourf(lulu,cmap=cmap_var2,norm=cnorm_var2,levels=clev_var2, extend='both') #extend='max'
        cbar2 = fig.colorbar(contf,ax=ax, location='left', pad=-0.06, extend='both') # shrink=0.7, ax=[axes[1],axes[3], axes[5]]
        cbar2.set_label(second_var)
        # cbar.set_label("SK38 / -")
        meta_x = 0.24
    else:
        meta_x = 0.17

    if var=='grain_type':
        lulu = np.zeros((n_bar,n_bar))
        for nn,k in enumerate(col_dict.keys()):
            lulu[nn, :] = np.nan # k
        norm_bins = np.sort([*col_dict.keys()]) + 0.5
        norm_bins = np.insert(norm_bins, 0, np.min(norm_bins) - 1.0)

        norm = BoundaryNorm(norm_bins, n_bar, clip=True)
        fmt = FuncFormatter(lambda x, pos: LABELS_GRAIN_TYPE_BAR[::-1][norm(x)])
        diff = norm_bins[1:] - norm_bins[:-1]
        tickz = norm_bins[:-1] + diff / 2

        # contf = ax.contourf(lulu,cmap=cmap,norm=norm,levels=norm_bins) # just for colorbar
        contf = ax.contourf(lulu,cmap=cmap,norm=norm,levels=norm_bins,hatches=HATCHES_GRAIN_TYPE_BAR[::-1], alpha=var_alpha) # just for colorbar
        cbar = fig.colorbar(contf, ax=ax, format=fmt, ticks=tickz,location='left', pad=0.01) # shrink=0.7, ax=[axes[1],axes[3], axes[5]]
        cbar.ax.grid(visible=False)
    else:
        n_var = 9
        lulu = np.zeros((n_var,n_var))
        for nn in range(0,n_var):
            lulu[nn, :] = np.nan # nn
        contf = ax.contourf(lulu,cmap=cmap_var,norm=cnorm_var,levels=clev_var, extend='both') #extend='max'
        cbar = fig.colorbar(contf,ax=ax, location='left', pad=0.01, extend='both') # shrink=0.7, ax=[axes[1],axes[3], axes[5]]
        cbar.set_label(var)
        # cbar.set_label("SK38 / -")
    
    """AXES""" 
    if DATE_RANGE[0] == 'NONE':
        ax.set_xlim(df_pro_list[0].date[0],df_pro_list[-1].date[0])
    else:
        # DATETIME_FORMAT  = '%Y-%m-%d' # +01:00
        # date_of_prof     = datetime.strptime(PROF_META['datetime'][0:10], DATETIME_FORMAT)
        ax.set_xlim(DATE_RANGE[0],DATE_RANGE[1])
    ax.set_ylim(0,np.max(h_max)+0.1)

    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_ylabel("height / m")

    # ax.xaxis.set_major_locator(###)
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    """Include Meta data in top left corner and save figure"""
    if DATETIME_STR!=None:
        meta_x = 0.13
        meta_y = 0.89
    else:
        meta_y = 0.94
    header_str = 'Location:      ' + meta_dict['StationName'] + '\nElevation:     ' + meta_dict['Altitude'] + \
                'm\nSlope Angle: ' + str(int(float(meta_dict['SlopeAngle']))) + '°\nAspect:         ' + str(int(float(meta_dict['SlopeAzi'])))  + '°'
    fig.text(meta_x,meta_y,header_str,horizontalalignment='left',
             verticalalignment='top', fontsize=10) # ma='left'
    
    # --- Save figure --- #
    if DATETIME_STR==None:
        fig_title = 'snp-evo-' + meta_dict['StationName'] + '.png'
    else:
        fig_title = 'snp-evo-and-profile-' + meta_dict['StationName'] + '.png'
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir,fig_title), facecolor='w', edgecolor='w',
                format='png', dpi=150)

    end_plotting = time.time()
    print('Visualization of snowpack evolution completed in {}s'.format(int(end_plotting-end_readin)))


def snowpro(config_file,pro_file=None):
    """A SNOWPACK output (.pro file) visualization tool.
    
    Arguments:
        config_file (str):  Path to configuration (ini) file   
        pro_file (str):     Path to PRO file
    """

    config = configparser.ConfigParser()
    config.read(config_file)

    os.makedirs(config.get('SNOWPRO','OUTPUT_DIR'))

    if pro_file != None:
        config['SNOWPRO','PRO_FILE_PATH'] = pro_file

    if config.get('SNOWPRO','PLOT_SNP_EVO') =='TRUE':
        DATE_RANGE = [config.get('SNOWPRO-EVO', 'START_DATE'), config.get('SNOWPRO-EVO','END_DATE')]
        plot_snp_evo(config.get('SNOWPRO','PRO_FILE_PATH'), output_dir=config.get('SNOWPRO','OUTPUT_DIR'), var=config.get('SNOWPRO-EVO','VAR'), res=config.get('SNOWPRO-EVO','RESOLUTION'), 
                        second_var=config.get('SNOWPRO-EVO','SECOND_VAR'), COLOR_SCHEME=config.get('SNOWPRO','COLOR_SCHEME'), DATE_RANGE=DATE_RANGE)
    
    if config.get('SNOWPRO','PLOT_PROFILE')=='TRUE':
        DATETIME = config.get('SNOWPRO-PROF', 'DATETIME')
        plot_single_profile(config.get('SNOWPRO','PRO_FILE_PATH'), config.get('PROFILE', 'DATETIME'))

    if config.get('SNOWPRO', 'PLOT_SNP_EVO_AND_PROFILE')=='TRUE':
        DATE_RANGE = [config.get('SNOWPRO-EVO', 'START_DATE'), config.get('SNOWPRO-EVO', 'END_DATE')]
        plot_snp_evo(config.get('SNOWPRO','PRO_FILE_PATH'), output_dir=config.get('SNOWPRO','OUTPUT_DIR'), DATETIME_STR=config.get('PROFILE','DATETIME'), var=config.get('SNOWPRO-EVO', 'VAR'), res=config.get('SNOWPRO-EVO', 'RESOLUTION'), 
                        second_var=config.get('SNOWPRO-EVO','SECOND_VAR'), COLOR_SCHEME=config.get('SNOWPRO','COLOR_SCHEME'), DATE_RANGE=DATE_RANGE)
        # plot_snp_evo(config.get('SNOWPRO','PRO_FILE'), DATETIME_STR=config.get('PROFILE','DATETIME'), var=config.get('SNOWPRO-evo', 'VAR'), res=config.get('SNOWPRO-evo', 'RESOLUTION'), 
        #                 second_var='NONE', COLOR_SCHEME=config.get('SNOWPRO','COLOR_SCHEME'), DATE_RANGE=DATE_RANGE)


if __name__ == "__main__":
    """
    import debugpy
    debugpy.listen(5678)
    print('Waiting for debugger!')
    debugpy.wait_for_client()
    print('Attached!')
    """

    if os.path.exists('latex_template.mplstyle'):
        plt.style.use('latex_template.mplstyle')

    args = sys.argv[1:]
    if (os.path.isfile(args[0])):
        if len(args)>1: 
            snowpro(args[0],args[1])
        else:
            snowpro(args[0])
    else:
        name = str(args[0])
        print('File ({}) not found.'.format(name))