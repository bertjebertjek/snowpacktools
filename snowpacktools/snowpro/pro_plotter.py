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
from snowpacktools.snowpro import snowpro, pro_helper, instability_rfm_mayer

latex_template_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "latex_template.mplstyle")
if os.path.exists(latex_template_path):
    plt.style.use(latex_template_path)

    
def snp_evo(pro_path: str, out_path: str=None, DATETIME_STR=None, var: str="grain_type", second_var: str=None,
            res: str="1h", color_scheme: str="SARP", height_max: float=None, date_range=None):
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
        var (str, optional): 
            Primary variable to visualize (default: "grain_type").
        second_var (str, optional): 
            Secondary variable to overlay (e.g., "Sk38", "Sn38", "Punstable", "ccl"). Defaults to None.
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
    if var=='grain_type':
        col_dict_labels     = dict(zip(LABELS_GRAIN_TYPE, COLORS_GRAIN_TYPE))
        hatches_dict_labels = dict(zip(LABELS_GRAIN_TYPE, HATCHES_GRAIN_TYPE))
    else:
        if var == 'Punstable':
            profs     = instability_rfm_mayer.calc_punstable(profs)
            cmap_var  = pro_helper.get_Punstable_cmap()
            var_ticks = np.arange(0,1.1,0.1)
        elif var == "ccl":
            cmap_var  = pro_helper.get_ccl_cmap()
            var_ticks = np.arange(0,1.1,0.1)
        elif var == 'Sk38' or var == 'Sn38':
            cmap_var  = pro_helper.get_sk38_cmap()
            var_ticks = np.arange(0,1.6,0.1)
        else:
            cmap_var = plt.get_cmap('plasma_r')
        clev_var = np.linspace(RANGE_DICT[var][0],RANGE_DICT[var][1],100) # 11
        # cnorm_var = BoundaryNorm(boundaries=clev_var, ncolors=cmap_var.N, clip=True)

    plot_second_var  = False
    var_alpha        = 1
    second_vars      = ['Sk38','Sn38','Punstable','ccl']
    hatch_second_var = ''
    if second_var in second_vars:
        plot_second_var = True
        var_alpha=0.33
        if second_var == 'Sk38' or second_var == 'Sn38':
            cmap_var2        = pro_helper.get_sk38_cmap()
            second_var_ticks = np.arange(0,1.6,0.1)
        elif second_var == 'Punstable':
            profs            = instability_rfm_mayer.calc_punstable(profs)
            cmap_var2        = pro_helper.get_Punstable_cmap()
            second_var_ticks = np.arange(0,1.1,0.1)
        else:
            #cmap_var2  = plt.get_cmap('gist_gray')
            cmap_var2 = pro_helper.get_whiteout_cmap(reverse=True)
            
        clev_var2 = np.linspace(RANGE_DICT[second_var][0],RANGE_DICT[second_var][1],100) # 11 discrete colorbar
        # cnorm_var2 = BoundaryNorm(boundaries=clev_var2, ncolors=cmap_var2.N, clip=False)

    if DATETIME_STR is None:
        # fig, (ax_cbar,ax) = plt.subplots(1,2,figsize=(10,5),gridspec_kw={'width_ratios': [0.5,11]})
        fig, ax = plt.subplots(1,1,figsize=(9.5,5))
    else:
        fig, (ax,ax_prof) = plt.subplots(1,2,figsize=(11,5),sharey=True,gridspec_kw={'width_ratios': [11,4]})
    # ax_cbar.axis('off')

    h_max = []
    dates = []
    for i, ts in enumerate(profs):
        prof = profs[ts]
        dates.append(ts)
        if len(prof['height']>0):
            h_max.append(prof['height'][-1])
                
            if var=='grain_type':
                cols   =[]
                hatches=[]
                for row in range(0,len(prof['graintype'])):
                    cols.append(col_dict_labels[prof['graintype'][row][0]])
                    hatches.append(hatches_dict_labels[prof['graintype'][row][0]])
                ax.bar(ts, prof['thickness'], width=w, bottom=prof['bottom'], align='edge', color=cols, hatch=hatches, alpha=var_alpha) # label=labels[i])
            else:
                thickness = np.where(prof['RTA']>=rta, prof['thickness'], np.nan)
                bottom    = np.where(prof['RTA']>=rta, prof['bottom'],    np.nan)

                var_data = (prof[var]-RANGE_DICT[var][0])/(RANGE_DICT[var][1]-RANGE_DICT[var][0])
                cols     = cmap_var(var_data)
                ax.bar(ts, prof['height'][-1], width=w, bottom=0, align='edge', color='lightgrey', alpha=1)
                ax.bar(ts, thickness, width=w, bottom=bottom, align='edge', color=cols, alpha=1)
            
            if plot_second_var:
                """Filter data with RTA"""
                thickness = np.where(prof['RTA']>=rta, prof['thickness'], np.nan)
                bottom    = np.where(prof['RTA']>=rta, prof['bottom'],    np.nan)
                hatches   = np.where(prof['RTA']>=rta, hatch_second_var, np.nan)

                var_data2 = np.where(prof['RTA']>=rta, prof[second_var],  np.nan)
                var_data2 = (var_data2-RANGE_DICT[second_var][0])/(RANGE_DICT[second_var][1]-RANGE_DICT[second_var][0])
                cols2 = cmap_var2(var_data2)
                ax.bar(ts, thickness, width=w, bottom=bottom, align='edge', color=cols2, hatch=hatches, alpha=1)
        else:
            h_max.append(0)

    if height_max is not None:
        height_max = 1.1*np.max(h_max)
    ax.set_ylim(0,height_max)

    ## Line along snow surface
    ax.plot(dates,h_max,ds='steps-post',lw=0.8,color='black', ls='--', alpha=0.67)

    ## Hardness profile to the right
    if DATETIME_STR is not None:
        ax_prof, DATETIME_STR = single_profile(pro_path, DATETIME_STR, height_max=height_max, color_scheme=color_scheme, ax=ax_prof)
        datetime_format = '%Y-%m-%dT%H:%M'
        time_of_profile = datetime.strptime(DATETIME_STR, datetime_format)
        ax.axvline(x=time_of_profile,ymin=-0.1, ymax=1.1, color='black', lw=1.5, ls='--')

    ## COLORBAR (Norm, bins, formatter, ticks - lots of stuff to make colorbar look nice)
    if second_var is not None:
         # - Colorbar for second layer - # 
        n_var = 9
        lulu = np.zeros((n_var,n_var))
        for nn in range(0,n_var):
            lulu[nn, :] = np.nan # nn
        contf = ax.contourf(lulu, cmap=cmap_var2, levels=clev_var2, hatches=hatch_second_var)
        # cbar2 = fig.colorbar(contf,ax=ax_cbar, location='left', ticks=second_var_ticks, pad=0, fraction=1, aspect=30, extend='both') # shrink=0.7, ax=[axes[1],axes[3], axes[5]]
        cbar2 = fig.colorbar(contf,ax=ax, location='left', ticks=second_var_ticks, pad=0.02, extend='both') # shrink=0.7
        cbar2.set_label(second_var)

        """Modify grain type legend for SARPGR"""
        LABELS_GRAIN_TYPE  = ['PP(gp), DF','SH, DH','FC(xr), RG','MF(cr), IF']
        COLORS_GRAIN_TYPE  = ['#ffde00','#95258f','#dacef4','#d5ebb5']
        HATCHES_GRAIN_TYPE = ['','','','']
        pro_helper.add_custom_legend(ax, LABELS_GRAIN_TYPE, COLORS_GRAIN_TYPE, HATCHES_GRAIN_TYPE, x=0.188, y=1.03, width=0.028, height=0.025, spacing=0.17, alpha=var_alpha)
    else:
        if var=='grain_type':
            pro_helper.add_custom_legend(ax, LABELS_GRAIN_TYPE[1:], COLORS_GRAIN_TYPE[1:], HATCHES_GRAIN_TYPE[1:], x=0.015, y=1.03, width=0.028, height=0.025, spacing=0.092, alpha=var_alpha)
        else:
            n_var = 9
            lulu = np.zeros((n_var,n_var))
            for nn in range(0,n_var):
                lulu[nn, :] = np.nan # nn
            contf = ax.contourf(lulu,cmap=cmap_var,levels=clev_var) #extend='max'
            # cbar = fig.colorbar(contf,ax=ax_cbar, location='left', ticks=var_ticks, fraction=1, extend='both')
            cbar = fig.colorbar(contf,ax=ax, location='left', ticks=var_ticks, pad=0.02, extend='both')
            cbar.set_label(var)
    
    ## Axes and labels
    if date_range is None:
        ax.set_xlim(dates[0],dates[-1])
    else:
        DATETIME_FORMAT  = '%Y-%m-%d' # +01:00
        d0 = datetime.strptime(date_range[0], DATETIME_FORMAT)
        d1 = datetime.strptime(date_range[1], DATETIME_FORMAT)
        ax.set_xlim(d0,d1)

    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.yaxis.set_tick_params(labelright=True)
    ax.set_ylabel("height / cm")
    
    date_fmt = DateFormatter("%b-%d")
    ax.xaxis.set_major_formatter(date_fmt)
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    ## Include Meta data in top left corner and save figure
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
        if DATETIME_STR is None:
            out_path = filename + '_snp_evo.png'
        else:
            out_path = filename + '_snp_evo_and_profile.png'

    fig.tight_layout()
    print(f"[i]  Saving figure {out_path}...")
    fig.savefig(out_path, facecolor='w', edgecolor='w', format='png', dpi=150, bbox_inches='tight')
    print(f"[i]  Visualization completed in {time.time()-start_time}s")


def single_profile(pro_path, DATETIME_STR, height_max=None, color_scheme="SARP", ax=None, out_path=None):
    """
    Plots a hardness profile of a snow profile from a PRO file for a specific timestamp including a temperature
    profile.

    Arguments:
        pro_path (str): 
            Path to the PRO file containing snowpack data.
        DATETIME_STR (str): 
            Timestamp to extract and plot the snow profile. Format: '%Y-%m-%dT%H:%M'.
        height_max (float, optional): 
            Maximum height of the plot in centimeters. If None, height is auto-scaled. Defaults to None.
        color_scheme (str, optional): 
            Color scheme for visualizing grain types. Defaults to "SARP".
        ax (matplotlib.axes.Axes, optional): 
            Existing axes to plot on. If None, a new figure and axes are created. Defaults to None.
        out_path (str, optional): 
            Path to save the output plot. If None, the figure will be saved in the current directory with 
            a default name.

    Returns (optional):
        matplotlib.axes.Axes: 
            The axes object (if `ax` is provided) and the DATETIME_STR string
    """


    profs, meta_dict = snowpro.read_pro(pro_path)

    LABELS_GRAIN_TYPE, COLORS_GRAIN_TYPE, HATCHES_GRAIN_TYPE, LABELS_GRAIN_TYPE_BAR, COLORS_GRAIN_TYPE_BAR, HATCHES_GRAIN_TYPE_BAR = pro_helper.get_grain_type_colors(color_scheme)

    ## Get closest timestamp in PRO file to DATETIME_STR
    datetime_format  = '%Y-%m-%dT%H:%M'    
    ts = datetime.strptime(DATETIME_STR, datetime_format)

    if ts in profs.keys():
        print(f"[i]  Timestamp {ts} found. Hardness profile will be plotted.")
    else:
        ts_query = ts
        ts = get_closest_ts(ts, sorted(profs.keys()))
        print(f"[i]  Timestamp {ts_query} not found instead using {ts}.")
        DATETIME_STR = datetime.strftime(ts, datetime_format)
    prof = profs[ts]

    hand_hardness_param_needed = True
    if np.max(abs(prof['hand hardness'])) > 10:
        hand_hardness_param_needed = False
    hand_hardness_dict, tickz_hh, tick_labels_hh = pro_helper.get_hand_hardness_N_dict()
    prof['hand_hardness_N'] = prof['hand hardness'] # Just initialisation
    ZERO_HH_VAL = 50

    if ax is None:
        fig, ax = plt.subplots(1,1,figsize=(4.5,5))
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
    hardness_last_step = None
    for ilayer in range(0,len(prof['graintype'])):
        cols.append(col_dict_labels[prof['graintype'][ilayer][0]])
        hatches.append(hatches_dict_labels[prof['graintype'][ilayer][0]])
        hardness = prof['hand hardness'][ilayer]
        if not hardness or np.isnan(hardness):
            hardness = hardness_last_step # error if no hardness at all
        else:
            hardness_last_step = hardness
        if hand_hardness_param_needed:
            prof['hand_hardness_N'][ilayer] = hand_hardness_dict[hardness]
    bar_plot   = ax.barh(prof['bottom'], prof['hand_hardness_N'], height=prof['thickness'], align='edge', color=cols, hatch=hatches) # label=labels[i])
    bar_plot_p = ax.barh(prof['bottom'], ZERO_HH_VAL, height=prof['thickness'], align='edge', color=cols, hatch=hatches)
        

    ## COLORBAR (Norm, bins, formatter, ticks - lots of stuff to make colorbar look nice)
    if fig is not None:
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

    ## Temperature axis
    ax_t = ax.twiny()
    ax_t.plot(prof['temperature'], prof['height'], color='#DC143C',lw=1)
    ax_t.grid(visible=False)
    ax_t.xaxis.tick_bottom()
    ax_t.xaxis.set_label_position("bottom")
    ax_t.set_xlabel("snow temperature / °C", color='#DC143C')
    ax_t.tick_params(axis='x', colors='#DC143C')
    ax_t.set_xlim(min(prof['temperature'])-2,0)

    ## AXES
    XLIM = [-1100,50] 
    ax.set_xlim(XLIM)
    if fig is not None:
        if height_max is not None:
            ax.set_ylim(0,float(height_max))
        else:
            ax.set_ylim(0,1.1*prof['height'][-1])
        ax.set_ylabel("height / cm")
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
    # tickz       = [-1000,-500,-250,-100,-20]
    # tick_labels = ['K','P','1F','4F','F']
    ax_hh.set_xticks(tickz_hh)
    ax_hh.set_xticklabels(tick_labels_hh)
    ax_hh.tick_params(axis="x",direction="in", pad=-18)

    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.set_xlabel("hand hardness / N")

    ax.yaxis.set_minor_locator(AutoMinorLocator())

    if fig is not None:
        meta_x = 0.04
        meta_y = 0.92
        header_str   = 'Location:' + '\nElevation:' + '\nSlope Angle:' + '\nAspect:'
        header_str_2 = meta_dict['StationName'] + '\n' + meta_dict['Altitude'] + 'm\n' + str(int(float(meta_dict['SlopeAngle']))) + '°\n' + str(int(float(meta_dict['SlopeAzi'])))  + '°'
        ax.text(meta_x,       meta_y, header_str,    horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, fontsize=10)
        ax.text(meta_x + 0.28, meta_y, header_str_2, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, fontsize=10)

        # --- Save figure --- #
        if out_path is None:
            out_path = pro_path.split("/")[-1].split(".")[0] + '_snow_profile.png'
        fig.tight_layout()
        print(f"[i]  Saving figure {out_path}")
        fig.savefig(out_path, facecolor='w', edgecolor='w', format='png', dpi=150)
    else:
        ax.text(0.04,0.93,DATETIME_STR,horizontalalignment='left',
                verticalalignment='top', fontsize=10, transform=ax.transAxes) # ma='left'
        return ax, DATETIME_STR
    

def get_closest_ts(ts, ts_list):
    """Get closest ts in a presorted list of timestamps:  ts_list = sorted(index)"""
    i = bisect_left(ts_list, ts)
    return min(ts_list[max(0, i-1): i+2], key=lambda t: abs(ts - t))