import os
from datetime import datetime
import numpy as np
# import xarray
import time
from bisect import bisect_left

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.ticker import AutoMinorLocator, FuncFormatter
from matplotlib.dates import DateFormatter
import matplotlib.dates as mdates
from snowpacktools.snowpro import snowpro, pro_helper, instability_rfm_mayer

latex_template_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "latex_template.mplstyle")
if os.path.exists(latex_template_path):
    plt.style.use(latex_template_path)

    
def plt_dry_instability(pro_path: str, out_path: str=None, DATETIME_STR=None, res: str="1h", 
                      color_scheme: str="SARP", height_max: float=None, date_range=None):
    """
    Plots snowpack evolution from a PRO file, visualizing different variables or grain types. Supports overlaying
    variables such as instability, grain type, or specific snow characteristics.

    Arguments:
        pro_path (str): 
            Path to the PRO file to be analyzed and visualized.
        out_path (str, optional): 
            Path to save the output plot. If None, the figure will be saved in the current directory.
        DATETIME_STR (str, optional): 
            Specific timestamp to overlay the hardness profile. Format: '%Y-%m-%dT%H:%M'. Defaults to None.
        res (str, optional): 
            Temporal resolution for the analysis (e.g., "1h"). Defaults to "1h".
        color_scheme (str, optional): 
            Color scheme for visualizations. Defaults to "SARP".
        height_max (float, optional): 
            Maximum height of the plot in centimeters. If None, height is auto-scaled. Defaults to None.
        date_range (list, optional): 
            Start and end dates for the plot (e.g., ['2023-12-01', '2023-12-31']). Defaults to None.
    """

    ## Read .pro file and get bar width
    start_time= time.time()
    profs, meta_dict = snowpro.read_pro(pro_path, res=res, keep_soil=False, consider_surface_hoar=True)

    w, __ = pro_helper.set_resolution(res)
    LABELS_GRAIN_TYPE, COLORS_GRAIN_TYPE, HATCHES_GRAIN_TYPE, _, _, _ = pro_helper.get_grain_type_colors(color_scheme)
    RANGE_DICT = pro_helper.get_range_dict()

    ## Color map and preprocessing
    rta = 0.8

    col_dict_labels     = dict(zip(LABELS_GRAIN_TYPE, COLORS_GRAIN_TYPE))
    hatches_dict_labels = dict(zip(LABELS_GRAIN_TYPE, HATCHES_GRAIN_TYPE))

    var = "Punstable"
    if var == "Punstable":
        profs     = instability_rfm_mayer.calc_punstable(profs)
        cmap_var  = pro_helper.get_Punstable_cmap()
        var_ticks = np.arange(0,1.1,0.1)
    elif var == "Sk38":
        cmap_var  = pro_helper.get_sk38_cmap()
        var_ticks = np.arange(0,1.6,0.1)
    clev_var = np.linspace(RANGE_DICT[var][0],RANGE_DICT[var][1],100) # 11

    fig, axes = plt.subplots(3,2,figsize=(11,12), sharex=True, gridspec_kw={'hspace': 0.03, 'wspace': 0.03, 'height_ratios': [4,4,3], 'width_ratios':[1,10]})
    
    axes[0,0].axis('off')
    axes[1,0].axis('off')
    axes[2,0].axis('off')

    ax = axes[0,1]
    ax_pu = axes[1,1]
    ax_pucb = axes[1,0]
    ax_insta = axes[2,1]
    # ax_sk = axes[]

    dates = mdates.date2num(sorted(profs.keys()))
    ax_insta.set_ylim(0,1)
    ax_insta.set_ylabel("P$_{unstable,max}$ (RTA > 0.8)")
    ax_insta.grid(False)
    ax_insta.yaxis.label.set_color('red')
    ax_insta.tick_params(axis='y', colors='red')

    ax_depth = ax_insta.twinx()
    ax_depth.set_ylabel("depth (weakest layer) / cm")

    line_y = [0,0.7,0.78]
    line_h = [0.7,0.08,0.22]
    line_c = ['#85CC6F','#FF9F00','#d1001f']
    ax_insta.barh(line_y, dates[-1]-dates[0], height=line_h, left=dates[0], color=line_c, align='edge', alpha=0.2)

    h_max = []
    dates = []
    inst_index_max = []
    inst_index_max_depth = []
    for ts in profs:
        prof = profs[ts]
        dates.append(ts)
        if len(prof['height']>0):
            h_max.append(prof['height'][-1])
                
            ## Grain type
            cols   =[]
            hatches=[]
            for row in range(0,len(prof['graintype'])):
                cols.append(col_dict_labels[prof['graintype'][row][0]])
                hatches.append(hatches_dict_labels[prof['graintype'][row][0]])
            ax.bar(ts, prof['thickness'], width=w, bottom=prof['bottom'], align='edge', color=cols, hatch=hatches, alpha=1) # label=labels[i])
            
            ## Instability index: Punstable, SK38, ...
            thickness = np.where(prof['RTA']>=rta, prof['thickness'], np.nan)
            bottom    = np.where(prof['RTA']>=rta, prof['bottom'],    np.nan)

            var_data = (prof[var]-RANGE_DICT[var][0])/(RANGE_DICT[var][1]-RANGE_DICT[var][0])
            cols     = cmap_var(var_data)
            ax_pu.bar(ts, prof['height'][-1], width=w, bottom=0, align='edge', color='lightgrey', alpha=1)
            ax_pu.bar(ts, thickness, width=w, bottom=bottom, align='edge', color=cols, alpha=1)
            
            ## Instability evolution
            inst_index = np.where(prof['RTA']>=rta, prof['Punstable'], np.nan)
            inst_index_max.append(np.nanmax(inst_index))
            inst_index_max_depth.append(prof['height'][-1] - prof['height'][np.argmax(inst_index)])

        else:
            h_max.append(0)

    ## Instability evolution
    ax_insta.scatter(dates, inst_index_max, s=7, marker='D', facecolors='red', edgecolors='none')
    ax_depth.scatter(dates, inst_index_max_depth, s=6, c='k', marker='x')
    
    ax_depth.set_ylim(0,120)

    if height_max is not None:
        height_max = 1.1*np.max(h_max)
    ax.set_ylim(0,height_max)
    ax_pu.set_ylim(0,height_max)

    ## Line along snow surface
    ax.plot(dates,h_max,ds='steps-post',lw=0.8,color='black', ls='--', alpha=0.67)
    ax_pu.plot(dates,h_max,ds='steps-post',lw=0.8,color='black', ls='--', alpha=0.67)


    ## Vertical dashed line for end of nowcast
    if DATETIME_STR is not None:
        datetime_format = '%Y-%m-%dT%H:%M'
        time_of_profile = datetime.strptime(DATETIME_STR, datetime_format)
        ax.axvline(x=time_of_profile,ymin=-0.1, ymax=1.1, color='black', lw=1.5, ls='--')
        ax_pu.axvline(x=time_of_profile,ymin=-0.1, ymax=1.1, color='black', lw=1.5, ls='--')
        ax_insta.axvline(x=time_of_profile,ymin=0, ymax=1, color='black', lw=1.5, ls='--')


    ## Colorbar grain type (Norm, bins, formatter, ticks - lots of stuff to make colorbar look nice)
    pro_helper.add_custom_legend(ax, LABELS_GRAIN_TYPE[1:], COLORS_GRAIN_TYPE[1:], HATCHES_GRAIN_TYPE[1:], x=0.015, y=1.03, width=0.028, height=0.025, spacing=0.092, alpha=1)
    
    ## Colorbar instability
    contf = ax_pucb.contourf(np.empty((9,9)),cmap=cmap_var,levels=clev_var) #extend='max'
    cbar = fig.colorbar(contf,ax=ax_pucb, location='left', ticks=var_ticks,fraction=1, extend='both')
    cbar.set_label("P$_{unstable}$ (RTA > 0.8)")
    
    ## Axes and labels
    if date_range is None:
        ax.set_xlim(dates[0],dates[-1])
    else:
        DATETIME_FORMAT  = '%Y-%m-%d' # +01:00
        d0 = datetime.strptime(date_range[0], DATETIME_FORMAT)
        d1 = datetime.strptime(date_range[1], DATETIME_FORMAT)
        ax.set_xlim(d0,d1)

    for ax0 in axes[0:2,1]:
        ax0.yaxis.tick_right()
        ax0.yaxis.set_label_position("right")
        ax0.yaxis.set_tick_params(labelright=True)
        ax0.set_ylabel("height / cm")
    
        date_fmt = DateFormatter("%b-%d")
        ax0.xaxis.set_major_formatter(date_fmt)
        ax0.yaxis.set_minor_locator(AutoMinorLocator())

    ## Include Meta data in top left corner of stratigraphy
    meta_x = 0.015
    meta_y = 0.975
    header_str   = 'Location:' + '\nElevation:' + '\nSlope Angle:' + '\nAspect:'
    header_str_2 = meta_dict['StationName'] + '\n' + meta_dict['Altitude'] + 'm\n' + str(int(float(meta_dict['SlopeAngle']))) + '°\n' + str(int(float(meta_dict['SlopeAzi'])))  + '°'
    ## cmb-prezis ##
    ax.text(meta_x,        meta_y, header_str,   horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, fontsize=10) # ma='left'
    ax.text(meta_x + 0.12, meta_y, header_str_2, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, fontsize=10) # ma='left'
    
    ## Save figure
    if out_path is None:
        filename = pro_path.split("/")[-1].split(".")[0]
        out_path = filename + '_dry_instability.png'

    fig.tight_layout()
    print(f"[i]  Saving figure {out_path}...")
    fig.savefig(out_path, facecolor='w', edgecolor='w', format='png', dpi=150, bbox_inches='tight')
    print(f"[i]  Visualization completed in {time.time()-start_time}s")