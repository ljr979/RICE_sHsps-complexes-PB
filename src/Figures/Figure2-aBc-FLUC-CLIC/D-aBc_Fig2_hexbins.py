"""#plots the hexbin plots for figure 2. This is a hexbin plot for all complexes at EACH timepoint, with the y axis being # of subunits of hsp, and the x axis being # of subunits of the matching client molecule.  (For aB-c)
Returns:
    hexbinplot
"""

import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
from loguru import logger
import random 
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

def read_counts(input_folder):
    """reads in the df which is formatted with matched molecule counts, pivoted to plot. concatinates multiple if necessary

    Args:
        input_folder (str): directory path within the repo where these counts live for each specific pair.

    Returns:
        df: dataframe with concatinated counts dataframes
    """
    counts_files_list=[filename for filename in os.listdir(f'{input_folder}') if 'counts_for_plotting' in filename]
    counts_collated=[]
    for filename in counts_files_list:
        count=pd.read_csv(f'{input_folder}{filename}')
        counts_collated.append(count)

    counts_collated=pd.concat(counts_collated)
    counts_collated.drop([col for col in counts_collated.columns.tolist() if 'Unnamed: 0' in col],axis=1)
    return counts_collated

def calc_ratios(counts_collated, pair):
    """calculate the ratio of shsp to client within each individual complex, at each timepoint. also counts the total number of complexes and puts this into a dataframe.

    Args:
        counts_collated (df): all the counts (matched) together
        pair (str  ): the proteins we are observing

    Returns:
        df, dict, dict (test, dict ratios, count_num_complexes: df with ratios in it, dictionary with the timepoint and the ratio of client/sHsp in it. dict with the number of complexes at each timepoint.
    """
    #removes any '0' counts
    better=counts_collated[counts_collated['client']>=1]
    better=better[better['hsp']>=1]
    #count the number of complexes
    total_num_complexes=len(better)

    test=[]
    for row, df in better.iterrows():
        row
        testo=pd.DataFrame(df).T
        #for each complex, how many hsp molecules and how many client molecules are there
        if testo['hsp'].values.tolist()[0] >= 1:
            testo['ratio']=testo['hsp']/testo['client']
            test.append(testo)
    test=pd.concat(test)
    #save
    test.to_csv(f'{output_folder}{pair}_ratios_added.csv')

    #this find the number of complexes (total) for each timepoint, and the median ratio of sHsp : client
    count_num_complexes={}
    dict_ratios={}
    for time, df in test.groupby('timepoint (min)'):
        med=df['ratio'].median()
        num_complexes_tp=len(df)
        dict_ratios[time]=med
        count_num_complexes[time]=num_complexes_tp
    return test, dict_ratios, count_num_complexes

def create_custom_colourmap(cmap, n_colors):
    """make a custom colourmap

    Args:
        cmap (str): the colourmap you want to adjust
        n_colors (int): number of shades you want in it

    Returns:
        colormap: colourmap with new range
    """
    new_colors = ['#ffffff'] + list(sns.color_palette(cmap, n_colors).as_hex())[10:40]

    #this turns it back into code hex for colourmap plotting
    
    cm = ListedColormap(new_colors)
   
    return cm

def plot_hexbin_kde(fig, axes, counts_collated, vmin, vmax, gridsize, extent, xlimo, ylimo, cm, cmap_axes):      
    """PLOTS the hexbins as however many timepoints as you had.

    Args:
        fig (fig): initialised out of function
        axes (axes): see fig
        counts_collated (df): dataframe with matched collated counts
        vmin (int ): this and vmax normalises the colour scale across all times so the intensity is relative between them to account for different number of observations
        vmax (int): max for normalisation
        xlimo (int): x limit
        ylimo (int): y limit
        cm (str): colourmap of choice
        cmap_axes (int): where to put the colormap (this has to be the bottom of the figure in this setup, so it is the number of rows in the figure, indexed to -1)

    Returns:
        fig, axes: the figure!
    """
    #now make it into subplots with normalised NUMBER OF OBSERVATIONS in each hex!!!!!
    # Setup figure with added gridspec to carve out an individual axes for the shared cbar
    for x, (timepoint, df) in enumerate(counts_collated.groupby('timepoint (min)')):
        x
        timepoint
        df
        # Add hexbin with set vmin, vax
        hexplot = axes[x].hexbin(df['client'], df['hsp'],gridsize=gridsize,
        extent=extent, cmap=cm, vmin=vmin, vmax=vmax)
        sns.kdeplot(df['client'], df['hsp'], color='darkgrey', linestyles='--', levels=np.arange(0, 1, 0.1), ax=axes[x])
        axes[x].set_xlim(0 , xlimo)
        axes[x].set_ylim(0, ylimo)
        axes[x].set(xlabel=None)
        axes[x].set_title(f'{timepoint} min', y=0.9, pad=-14)
        
    # add a shared colorbar in the last axes
    cb = fig.colorbar(hexplot, cax=axes[cmap_axes], orientation='horizontal')
    cb.set_label('Count')
    return fig, axes

def functional(calc_ratios, create_custom_colourmap, plot_hexbin_kde, output_folder, pair, cmap, ylimo, xlimo, gridsize, nrows, ncols, n_colors, vmin, vmax, cmap_axes, counts_collated):
    """runs the functions and outputs the hexbinplots for each pair of client/sHsp

    Args:
        calc_ratios (function): calculates the sHsp : client ratios
        create_custom_colourmap (function): creates a colourmap with limits you define
        plot_hexbin_kde (function): does the plotting
        output_folder (str): where to save the plots
        pair (str): the client and sHsp you're plotting
        cmap (str): the colours to make the colormap from
        ylimo (int): plotting limit (change according to the data you're plotting)
        xlimo (int): plotting limit (change according to the data you're plotting)
        nrows (int): number of timepoints in that dataset
        ncols (int): nrows -1
        n_colors (int): number of shades in your colourmap
        vmin (int): this is how you keep the same colour scale across timepoints
        vmax (int): this is how you keep the same colour scale across timepoints
        cmap_axes (int): _description_
        counts_collated (df): all the collated counts
    """
    counts_collated, dicto, complexes_count = calc_ratios(counts_collated, pair)
    #exclude last timepoint because FLUC does not have this same timepoint
    counts_collated = counts_collated[counts_collated['timepoint (min)']!=420.0]
    extent = [0, np.max(counts_collated['client']), 0, np.max(counts_collated['client'])]
    cm = create_custom_colourmap(cmap, n_colors)
    fig, axes = plt.subplots(nrows, ncols, figsize=(2,10), gridspec_kw={'height_ratios': [7, 7, 7, 7, 7, 0.5] })
    fig, axes = plot_hexbin_kde(fig, axes, counts_collated, vmin, vmax, gridsize, extent, xlimo, ylimo, cm, cmap_axes)
    #fig.savefig(f'{output_folder}{pair}.svg')
    #fig.savefig(f'{output_folder}{pair}.png')
    plt.show()


if __name__ == "__main__":
    
    #set font parameters
    font = {
        'family' : 'arial',
        'weight' : 'normal',
        'size'   : 12
    }
    plt.rc('font', **font)
    plt.rcParams['svg.fonttype'] = 'none'

    input_folder = 'data/Figures/Figure_2/D-matched_hexbins/aB-c/CLIC_aB-c/'
    output_folder = 'data/Figures/Figure_2/D-matched_hexbins/aB-c/'

    if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    counts_collated = read_counts(input_folder=input_folder)
    #specify the function inputs depending on the pair of proteins (colour, x and y limits, etc.)
    functional(calc_ratios, create_custom_colourmap, plot_hexbin_kde, output_folder, pair='CLIC and aB-c', cmap='Reds', ylimo=20, xlimo=20, gridsize=(30,20), nrows=6, ncols=1, n_colors=50, vmin=0, vmax=5, cmap_axes=5, counts_collated=counts_collated)

    #now for FLUC!    
    input_folder = 'data/Figures/Figure_2/D-matched_hexbins/aB-c/FLUC_aB-c/'
    output_folder = 'data/Figures/Figure_2/D-matched_hexbins/aB-c/'

    counts_collated = read_counts(input_folder=input_folder)
    functional(calc_ratios, create_custom_colourmap, plot_hexbin_kde, output_folder, pair='FLUC and aB-c', cmap='Purples', ylimo=40, xlimo=30, gridsize=(22,12), nrows=6, ncols=1, n_colors=50, vmin=0, vmax=10, cmap_axes=5, counts_collated=counts_collated)

