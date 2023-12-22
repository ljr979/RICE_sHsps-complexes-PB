import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
from loguru import logger
import random 
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
#plots the hexbin plots for figure 2. This is a hexbin plot for all complexes at EACH timepoint, with the y axis being # of subunits of hsp, and the x axis being # of subunits of the matching client molecule. 

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
    test.to_csv(f'{output_folder}{pair}stoichiometries_ratios_added.csv')

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

def plot_hexbin_kde(fig, axes, counts_collated, vmin, vmax, xlimo, ylimo, cm, cmap_axes):      
    #now make it into subplots with normalised NUMBER OF OBSERVATIONS in each hex!!!!!
    # Setup figure with added gridspec to carve out an individual axes for the shared cbar
    for x, (timepoint, df) in enumerate(counts_collated.groupby('timepoint (min)')):
        x
        timepoint
        df
        # Add hexbin with set vmin, vax
        hexplot = axes[x].hexbin(df['client'], df['hsp'], gridsize=30, cmap=cm, vmin=vmin, vmax=vmax)
        sns.kdeplot(df['client'], df['hsp'], color='darkgrey', linestyles='--', levels=np.arange(0, 1, 0.1), ax=axes[x])
        axes[x].set_xlim(0 , xlimo)
        axes[x].set_ylim(0, ylimo)
        axes[x].set(xlabel=None)
        axes[x].set_title(f'{timepoint} min', y=0.9, pad=-14)
        
    # add a shared colorbar in the last axes
    cb = fig.colorbar(hexplot, cax=axes[cmap_axes], orientation='horizontal')
    cb.set_label('Count')
    return fig, axes


if __name__ == "__main__":
        
    input_folder = 'data/Figures/Figure_2/D-matched_hexbins/aB-c/CLIC_aB-c/'
    output_folder = 'data/Figures/Figure_2/D-matched_hexbins/aB-c/'

    if not os.path.exists(output_folder):
            os.makedirs(output_folder)


    #set font parameters
    font = {
        'family' : 'arial',
        'weight' : 'normal',
        'size'   : 8
    }
    plt.rc('font', **font)
    plt.rcParams['svg.fonttype']= 'none'
    #change this for each pair of client:sHsp
    pair='CLIC and aB-c'
    #different colour for each pair
    cmap='Reds'
    #specify according to the data
    ylimo=20
    xlimo=20
    nrows=6
    ncols=1
    n_colors=50
    vmin=0
    vmax=5
    
    cmap_axes=nrows-1

    counts_collated=read_counts(input_folder=input_folder)

    counts_collated, dict_ratios, CLIC_complexes_count=calc_ratios(counts_collated, pair)
    #exclude last timepoint because FLUC does not have this same timepoint
    counts_collated=counts_collated[counts_collated['timepoint (min)']!=420.0]
    cm=create_custom_colourmap(cmap, n_colors)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5, 28), gridspec_kw={'height_ratios': [4, 4, 4, 4, 4, 0.5] })
    fig, axes=plot_hexbin_kde(fig, axes, counts_collated, vmin, vmax, xlimo, ylimo, cm, cmap_axes)
    fig.savefig(f'{output_folder}{pair}.svg')
    fig.savefig(f'{output_folder}{pair}.png')
    plt.show()
        



    #now for FLUC!    
    input_folder = 'data/Figures/Figure_2/D-matched_hexbins/aB-c/FLUC_aB-c/'
    output_folder = 'data/Figures/Figure_2/D-matched_hexbins/aB-c/'

    #change this for each pair of client:sHsp
    pair='FLUC and aB-c'

    cmap='Purples'

    ylimo=40
    xlimo=30
    nrows=6
    ncols=1
    n_colors=50
    vmin=0
    vmax=10
    cmap_axes=nrows-1

    counts_collated=read_counts(input_folder=input_folder)
    counts_collated, dict_ratios, FLUC_complexes_count=calc_ratios(counts_collated, pair)
    cm=create_custom_colourmap(cmap, n_colors)

    fig, axes = plt.subplots(nrows, ncols, figsize=(5, 28), gridspec_kw={'height_ratios': [4, 4, 4, 4, 4, 0.5], })
    fig, axes=plot_hexbin_kde(fig, axes, counts_collated, vmin, vmax, xlimo, ylimo, cm, cmap_axes)
    fig.savefig(f'{output_folder}{pair}.svg')
    fig.savefig(f'{output_folder}{pair}.png')
    plt.show()


