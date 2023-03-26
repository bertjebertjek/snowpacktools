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
from datetime import datetime

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.ticker import AutoMinorLocator, FuncFormatter
import matplotlib.image as image
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

from snowpacktools.snowpro import snowpro
from snowpacktools.snowpro import pro_helper


def plot_aps_and_profile_evolution(df_P, path_to_pro, output_path='output/', var='grain_type', res='1h', second_var='NONE', COLOR_SCHEME='IACS2',DATE_RANGE=['NONE','NONE']):
    """Visualization of avalanche problems and snowpack evolution in one figure with two axes.
    Plots snowpack evolution (PRO-file). Different variables or grain type can be visualized and overlayed.
    
    Arguments:
        path_to_pro (str):      Path to PRO file to be plotted
        output_dir (str):     
        DATETIME_STR (str):
        var: optional           Can be used for visualizing other variables than grain type
        res (str):
        second_var (str):
        COLOR_SCHEME (str):
        DATE_RANGE (list):
    Returns:
        Figure
    """

    ### Load mpl style
    this_dir, this_filename = os.path.split(__file__)
    latex_template_path = os.path.join(this_dir, "../snowpro/latex_template.mplstyle")
    if os.path.exists(latex_template_path):
        plt.style.use(latex_template_path)

    ### Load icons of APs
    icon_path       = "graphics/icons-aps-EAWS-colour/jpg"
    icon_newsnow    = os.path.join(this_dir, icon_path, 'Icon-Avalanche-Problem-New-Snow-EAWS.jpg')
    icon_newsnow    = image.imread(icon_newsnow)
    icon_drift      = os.path.join(this_dir, icon_path, 'Icon-Avalanche-Problem-Wind-Drifted-Snow-EAWS.jpg')
    icon_drift      = image.imread(icon_drift)
    icon_persistent = os.path.join(this_dir, icon_path, 'Icon-Avalanche-Problem-Persistent-Weak-Layer-EAWS.jpg')
    icon_persistent = image.imread(icon_persistent)
    icon_deep_pwl   = os.path.join(this_dir, icon_path, 'Icon-Avalanche-Problem-Deep-Persistent-Weak-Layer-EAWS.jpg')
    icon_deep_pwl   = image.imread(icon_deep_pwl)
    icon_wet        = os.path.join(this_dir, icon_path, 'Icon-Avalanche-Wet-Snow-EAWS.jpg')
    icon_wet        = image.imread(icon_wet)
    icon_gliding    = os.path.join(this_dir, icon_path, 'Icon-Avalanche-Problem-Gliding-Snow-EAWS.jpg')
    icon_gliding    = image.imread(icon_gliding)
    icon_list = [icon_newsnow, icon_drift, icon_persistent, icon_deep_pwl, icon_wet, icon_gliding]

    start_readin = time.time()
    df_pro_list_temp, meta_dict = snowpro.read_pro(path_to_pro,res=res)
    end_readin = time.time()
    print('[I]  Reading of PRO file completed in {}s'.format(int(end_readin-start_readin)))

    # Filter for certain resolution and time frame
    w, hours = pro_helper.set_resolution(res)
    LABELS_GRAIN_TYPE, COLORS_GRAIN_TYPE, HATCHES_GRAIN_TYPE, LABELS_GRAIN_TYPE_BAR, COLORS_GRAIN_TYPE_BAR, HATCHES_GRAIN_TYPE_BAR = pro_helper.get_grain_type_colors(COLOR_SCHEME)
    RANGE_DICT = pro_helper.get_range_dict()

    df_pro_list = []
    for df in df_pro_list_temp:
        if df.date.iloc[0].hour in hours:
        # if (df.iloc[0].date.strftime('%m.%d') > season_start) or (df.iloc[0].dates.strftime('%m.%d') < season_end):
            df_pro_list.append(df)

    # COLOR MAP AND PREPROCESSING
    if var=='grain_type':
        col_dict_labels     = dict(zip(LABELS_GRAIN_TYPE, COLORS_GRAIN_TYPE))
        hatches_dict_labels = dict(zip(LABELS_GRAIN_TYPE, HATCHES_GRAIN_TYPE))

        n_bar    = len(LABELS_GRAIN_TYPE_BAR)
        col_nums = np.arange(0,n_bar)
        col_dict = dict(zip(col_nums, COLORS_GRAIN_TYPE_BAR[::-1]))
        cmap     = ListedColormap([col_dict[x] for x in col_dict.keys()])
    
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

    # VISUALIZATION
    fig, ((ax0, ax),(ax1,ax_aps)) = plt.subplots(2,2,figsize=(12,7),sharex=True, gridspec_kw={'width_ratios':[1,15],'hspace':0.03,'wspace':0.03}) # 'height_ratios':[1,1]
    ax0.axis('off')
    ax1.axis('off')

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

    # Line along snow surface
    if second_var!='NONE':
        # h_max = np.where(h_max == np.nan, 0, h_max)
        ax.plot(dates,h_max,ds='steps-post',lw=0.8,color='black', ls='--', alpha=0.67)

    # COLORBAR (Norm, bins, formatter, ticks - lots of stuff to make colorbar look nice)
    if second_var!='NONE':
         # - Colorbar for second layer - # 
        n_var = 9
        lulu = np.zeros((n_var,n_var))
        for nn in range(0,n_var):
            lulu[nn, :] = np.nan # nn
        contf = ax.contourf(lulu,cmap=cmap_var2,norm=cnorm_var2,levels=clev_var2, extend='both') #extend='max'
        cbar2 = fig.colorbar(contf,ax=ax0, location='left', pad=-0.06, extend='both') # shrink=0.7, ax=[axes[1],axes[3], axes[5]]
        cbar2.set_label(second_var)
        # cbar.set_label("SK38 / -")

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
        cbar = fig.colorbar(contf, ax=ax0, format=fmt, ticks=tickz,location='left', fraction=1) # pad=0.01, shrink=0.7, ax=[axes[1],axes[3], axes[5]]
        cbar.ax.grid(visible=False)
    else:
        n_var = 9
        lulu = np.zeros((n_var,n_var))
        for nn in range(0,n_var):
            lulu[nn, :] = np.nan # nn
        contf = ax.contourf(lulu,cmap=cmap_var,norm=cnorm_var,levels=clev_var, extend='both') #extend='max'
        cbar = fig.colorbar(contf,ax=ax0, location='left', pad=0.01, extend='both') # shrink=0.7, ax=[axes[1],axes[3], axes[5]]
        cbar.set_label(var)
        # cbar.set_label("SK38 / -")
    
    # AXES 
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

    # Include Meta data in top left corner and save figure
    meta_x = 0.015
    meta_y = 0.98
    header_str = 'Location:      ' + meta_dict['StationName'] + '\nElevation:     ' + meta_dict['Altitude'] + \
                'm\nSlope Angle: ' + str(int(float(meta_dict['SlopeAngle']))) + '°\nAspect:         ' + str(int(float(meta_dict['SlopeAzi'])))  + '°'
    ax.text(meta_x,meta_y,header_str,horizontalalignment='left',
             verticalalignment='top', fontsize=10, transform=ax.transAxes) # ma='left'
    
    
    """Visusalize APs (second axis)"""
    df_P['napex_sele_natural'] = np.where(df_P['napex_sele_natural']==1, df_P['napex_sele_natural'], np.nan)
    df_P['papex_sele_natural'] = np.where(df_P['papex_sele_natural']==1, df_P['papex_sele_natural'], np.nan)
    df_P['dapex_sele_natural'] = np.where(df_P['dapex_sele_natural']==1, df_P['dapex_sele_natural'], np.nan)

    # LABELS_GRAIN_TYPE_BAR        = ['PP','DF','PPgp','SH','DH','FC(xr)','RG','MF','MFcr','IF']
    # COLORS_GRAIN_TYPE_BAR_IACS2  = ['#00FF00','#228B22','#696969','#FF00FF','#0000FF','#ADD8E6','#FFB6C1','#FF0000','#FF0000','#00FFFF']
    COLORS_APS = {'newSnow':'#00FF00','windSlab':'#228B22','deepPW':'#0000FF','PW':'#ADD8E6','wetSnow':'#FF0000'}
    LW_NATURAL = 1
    C_NATURAL  = "gold"
    ax_aps.set_ylim([0,6])

    ax_aps.bar(df_P['dy'], df_P['napex_sele_trigger'], width=0.75, bottom=5, color=COLORS_APS['newSnow'])
    ax_aps.bar(df_P['dy'], df_P['winex'], width=0.75, bottom=4, color=COLORS_APS['windSlab'])
    ax_aps.bar(df_P['dy'], df_P['papex_sele_trigger'], width=0.75, bottom=3, color=COLORS_APS['PW'])
    ax_aps.bar(df_P['dy'], df_P['dapex_sele_trigger'], width=0.75, bottom=2, color=COLORS_APS['deepPW'])
    ax_aps.bar(df_P['dy'], df_P['wapex_sele'], width=0.75, bottom=1, color=COLORS_APS['wetSnow'])

    ax_aps.bar(df_P['dy'], df_P['napex_sele_natural'], width=0.75, bottom=3, color=COLORS_APS['newSnow'], edgecolor = C_NATURAL, linewidth=LW_NATURAL)
    ax_aps.bar(df_P['dy'], df_P['papex_sele_natural'], width=0.75, bottom=3, color=COLORS_APS['PW'], edgecolor = C_NATURAL, linewidth=LW_NATURAL, label="Natural release")
    ax_aps.bar(df_P['dy'], df_P['dapex_sele_natural'], width=0.75, bottom=3, color=COLORS_APS['deepPW'], edgecolor = C_NATURAL, linewidth=LW_NATURAL)
    
    ax_aps.legend(loc='upper right')
    ### Add icons of avalanche problems
    icon_xpos = -0.03
    icon_ypos = np.linspace(0.085,0.915,6)[::-1]
    shade = 0.75
    pad = 0.08
    icon_zoom = 0.045
    
    for i,icon in enumerate(icon_list):
        imagebox = OffsetImage(icon, zoom=icon_zoom)
        ab = AnnotationBbox(imagebox, (icon_xpos,icon_ypos[i]), xycoords="axes fraction", frameon=True, pad=pad)
        ax_aps.add_artist(ab)
    ax_aps.yaxis.set_ticklabels([])

    # --- Save figure --- #    
    fig.tight_layout()
    print('[I]  Saving AP and snowpack evolution figure')
    fig.savefig(output_path, facecolor='w', edgecolor='w',
                format='png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    end_plotting = time.time()
    print('[I]  Visualization of APs and snowpack evolution completed in {}s'.format(int(end_plotting-end_readin)))