import os
from datetime import datetime
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
import re, sys

############################################
#
# Functions to read and plot SMET files
#
# Author: 2024, Bert Kruyt, NVE
############################################


def read_lines_after_keyword(filename, keyword, num_lines=None):
    """ Read text file line by line from """
    with open(filename, "r") as file:
        lines = file.readlines()

    # Find the keyword
    for i, line in enumerate(lines):
        if "field" in line:
            cols = line.split(" = ")[-1]
            # print(cols)
        if keyword in line:
            # Read lines after the keyword
            if num_lines:
                return (cols, lines[i+1:i+1+num_lines])  # Read a specific number of lines
            else:
                return (cols, lines[i+1:])  # Read all lines after the keyword

    return []  # Return an empty list if keyword is not found


def read_lines_between_keywords(file_path, start_keyword="[HEADER]", end_keyword="[DATA]"):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    start_index = None
    end_index = None
    delimiter2='\s+' # one or more spaces

    for i, line in enumerate(lines):
        if start_keyword in line:
            start_index = i
        if end_keyword in line:
            end_index = i
            break

    if start_index is not None and end_index is not None:
        data_lines = lines[start_index + 1:end_index]
        data_dict = {}
        for line in data_lines:
            key, value = line.strip().split('=', 1)
            # data_dict[key.strip()] = value.strip().split(delimiter2)
            data_dict[key.strip()] = re.split(r'\s+', value.strip() )
        return data_dict
    else:
        return None


### For meteio smet (input) smet files only
def SMET_in2df(filename, keyword='[DATA]', num_lines=None):
    """ 
    Read an input smet file and return a DataFrame 
    filename : str : input file name
    keyword : str : keyword to read below in the file
    num_lines : int : number of lines to read after the keyword (optional)

    
        return df : DataFrame : DataFrame of the data
    """
    
    
    # Read lines after the keyword    
    (col_names, lines_after_keyword) = read_lines_after_keyword(filename, keyword, num_lines=None) 

    delimiter1=" "
    col_names = col_names.strip().split(delimiter1)
    
    # Convert to DataFrame
    delimiter2="\t"

    data = [line.strip().split(delimiter2) for line in lines_after_keyword]

    # Create DataFrame with optional column names
    df = pd.DataFrame(data, columns=col_names if col_names else None)

    # Convert to float and datetime
    df[col_names[1:]] = df[col_names[1:]].astype(float)
    df['timestamp'] = pd.to_datetime(df['timestamp'] ) #, format='%Y-%m-%d %H:%M:%S')
    return df  


def SMET2df(filename, keyword='[DATA]', num_lines=None):
    """
    Read an input smet file and return a DataFrame and header dictionary

    Arguments:
    filename : str : input file name
    keyword : str : keyword to read below in the file
    num_lines : int : number of lines to read after the keyword (optional)

    
        return df , header : DataFrame : DataFrame of the data, header: dictionary of the header
    """
    
    
    # Read lines after the keyword    
    (col_names, lines_after_keyword) = read_lines_after_keyword(filename, keyword, num_lines=None) 

    delimiter1=" "
    col_names = col_names.strip().split(delimiter1)
    
    # Convert to DataFrame
    delimiter2="\t"
    data = [line.strip().split(delimiter2) for line in lines_after_keyword]
    if len(data[0]) == 1:
        delimiter2='\s+' # one or more spaces
        data = [line.strip().split(delimiter2) for line in lines_after_keyword]
        data = [re.split(r'\s+', line.strip()) for line in lines_after_keyword]
        if len(data[0]) == 1:
            delimiter2=' ' # one spaces
            data = [line.strip().split(delimiter2) for line in lines_after_keyword]

    # print(len(data[0]))
    
    if len(data[0]) > len(col_names):
        data = [line[:-1] for line in data]

        print(len(data[0]))


    # Create DataFrame with optional column names
    df = pd.DataFrame(data, columns=col_names if col_names else None)

    # Convert to float and datetime
    df[col_names[1:]] = df[col_names[1:]].astype(float)
    df['timestamp'] = pd.to_datetime(df['timestamp'] ) #, format='%Y-%m-%d %H:%M:%S')


    ### get the header
    header = read_lines_between_keywords(filename, start_keyword = '[HEADER]', end_keyword = keyword)


    return df  , header



#### Plot functions
def plot_vars(smet, vars_to_plot, start_date=None, end_date=None, header_dict=None, fig=None, ax=None, label=None, title=None):
    """ 
    Plot variables vars_to_plot from the SMET file 

    Arguments:
    smet : DataFrame : DataFrame of the SMET file with columns as variables
    vars_to_plot : list : list of variables to plot ( should be columns in the DataFrame). Note that when providing one var, it should still be a list: e.g. ['var']
    start_date : str : start date for the plot (default is None and will plot all dates)
    end_date : str : end date for the plot (default is None and will plot all dates)
    header_dict : dict : header dictionary with metadata ( from read_header function)
    fig : Figure : Figure object for plotting (optional)
    ax : Axes : Axes object for plotting (optional)
    label : str : label for the plot (optional)
    title : str : title for the plot (optional)
    """
    # fig, ax = plt.subplots( int(len(vars_to_plot)/2),2,figsize=(7*int(len(vars_to_plot)/2),4*2))
    if not fig:
        fig, ax = plt.subplots( int(np.ceil(len(vars_to_plot)/2)),2,figsize=(7*2, 4*np.ceil(len(vars_to_plot)/2)))

    if ( type(ax) == np.ndarray ) : #and ( len(ax.flatten()) >= len(vars_to_plot) ) :
        if  ( len(ax.flatten()) < len(vars_to_plot) ) :
            sys.exit('Provide more axes for the number of variables to plot')
        else:
            for v, var in enumerate(vars_to_plot):
                ax.flatten()[v].plot( smet['timestamp'], smet[var], label=label ) 
                ax.flatten()[v].set_title(var)
                ax.flatten()[v].grid(color='gray', linestyle='--', linewidth=0.5)
                if start_date and end_date:
                    ax.flatten()[v].set_xlim(pd.to_datetime(start_date), pd.to_datetime(end_date))
                # rotate xtixk labels
                ax.flatten()[v].set_xticks(ax.flatten()[v].get_xticks())
                ax.flatten()[v].set_xticklabels(ax.flatten()[v].get_xticklabels(), rotation=45, ha='right')
                # units from the header dictionary:        
                if header_dict:
                        ax.flatten()[v].set_ylabel( header_dict['plot_unit'][ smet.columns.get_loc(var) ]                )
                if title:
                    fig.suptitle(title)
                elif header_dict:
                    fig.suptitle(f"{header_dict['station_name'][0]} lat:{header_dict['latitude'][0]} lon:{header_dict['longitude'][0]} {header_dict['altitude'][0]}m")
                
                    
                # ...
                if label:
                    ax.flatten()[v].legend()
    elif len(vars_to_plot) > 1 and  ( type(ax) is not np.ndarray ) :
        sys.exit('Provide more axes for the number of variables to plot')
    else:  # i.e.  ( type(ax) is not np.ndarray ) so one plot, and one variable
        # for v, var in enumerate(vars_to_plot):
            var = vars_to_plot[0]
            ax.plot( smet['timestamp'], smet[var], label=label )
            ax.set_title(var)
            ax.grid(color='gray', linestyle='--', linewidth=0.5)
            if start_date and end_date:
                ax.set_xlim(pd.to_datetime(start_date), pd.to_datetime(end_date))
            # rotate xtixk labels
            ax.set_xticks(ax.get_xticks())
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
            # units from the header dictionary:        
            if header_dict:
                    ax.set_ylabel( header_dict['plot_unit'][ smet.columns.get_loc(var) ]           )
            if title:
                fig.suptitle(title)
            elif header_dict:
                fig.suptitle(f"{header_dict['station_name'][0]} lat:{header_dict['latitude'][0]} lon:{header_dict['longitude'][0]} {header_dict['altitude'][0]}m")
            
                
            # ...
            if label:
                ax.legend()

    plt.tight_layout()