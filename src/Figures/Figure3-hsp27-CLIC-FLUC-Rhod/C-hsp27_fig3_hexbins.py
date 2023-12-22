import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
from loguru import logger
import random 
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
#plotting figure 3 hexplots (at each timepoint, the number of subunits in each complex)

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
    better=counts_collated[counts_collated['client']>=1]
    better=better[better['hsp']>=1]
    total_num_complexes=len(better)
    test=[]
    for row, df in better.iterrows():
        row
        testo=pd.DataFrame(df).T
        if testo['hsp'].values.tolist()[0] >= 1:
            testo['ratio']=testo['hsp']/testo['client']
            test.append(testo)
    test=pd.concat(test)
    test.to_csv(f'{output_folder}{pair}stoichiometries_ratios_added.csv')

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


def plot_hexbin_kde(fig, axes, counts_collated, vmin, vmax, xlimo, ylimo, cm, cmap_axes, gridsize):      
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
    for x, (timepoint, df) in enumerate(counts_collated.groupby('timepoint (min)')):
        x
        timepoint
        df
        # Add hexbin with set vmin, vax
        hexplot = axes[x].hexbin(df['client'], df['hsp'], gridsize=gridsize, cmap=cm, vmin=vmin, vmax=vmax)
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
    font = {
        'family' : 'arial',
        'weight' : 'normal',
        'size'   : 16
    }
    plt.rc('font', **font)
    plt.rcParams['svg.fonttype']= 'none'

    input_folder = 'data/Figures/Figure_3/C-matched_hexbins/hsp27/Rhod_hsp27/'
    output_folder = 'data/Figures/Figure_3/C-matched_hexbins/hsp27/'

    if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    #change this for each pair of client:sHsp
    pair='Rhodanese and hsp27'

    cmap='Greens'
    #according to whatever you're comparing to, or what the max count is
    ylimo=25
    xlimo=30
    #the number of timepoints
    nrows=7
    #this makes it into a vertical plot, with 7 boxes in one column
    ncols=1

    #number of colours in your colormap (more shades)
    n_colors=50


    vmin=0
    vmax=30

    cmap_axes=nrows-1
    gridsize=20

    counts_collated=read_counts(input_folder=input_folder)
    counts_collated, dict_ratios, rhod_complexes_count=calc_ratios(counts_collated, pair)
    cm=create_custom_colourmap(cmap, n_colors)

    fig, axes = plt.subplots(nrows, ncols, figsize=(5, 35), gridspec_kw={'height_ratios': [4, 4, 4, 4, 4, 4, 0.5], })
    fig, axes=plot_hexbin_kde(fig, axes, counts_collated, vmin, vmax, xlimo, ylimo, cm, cmap_axes, gridsize)
    fig.savefig(f'{output_folder}{pair}.svg')
    plt.show()


    input_folder = 'data/Figures/Figure_3/C-matched_hexbins/hsp27/FLUC_hsp27/'
    output_folder = 'data/Figures/Figure_3/C-matched_hexbins/hsp27/'

    if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    #change this for each pair of client:sHsp
    pair='FLUC_hsp27'

    cmap='Purples'

    ylimo=40
    xlimo=30
    nrows=7
    ncols=1
    n_colors=50


    vmin=0
    vmax=30

    cmap_axes=nrows-1
    gridsize=20


    counts_collated=read_counts(input_folder=input_folder)
    counts_collated, dict_ratios, FLUC_complexes_count=calc_ratios(counts_collated, pair)
    cm=create_custom_colourmap(cmap, n_colors)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5, 35), gridspec_kw={'height_ratios': [4, 4, 4, 4, 4, 4, 0.5]})
    fig, axes=plot_hexbin_kde(fig, axes, counts_collated, vmin, vmax, xlimo, ylimo, cm, cmap_axes, gridsize)
    fig.savefig(f'{output_folder}{pair}.svg')
    plt.show()




    input_folder = 'data/Figures/Figure_3/C-matched_hexbins/hsp27/CLIC_hsp27/'
    output_folder = 'data/Figures/Figure_3/C-matched_hexbins/hsp27/CLIC_hsp27/'

    if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    #change this for each pair of client:sHsp
    pair='CLIC and hsp27'

    cmap='Reds'

    ylimo=10
    xlimo=16
    nrows=7
    ncols=1
    n_colors=50


    vmin=0
    vmax=8

    cmap_axes=nrows-1
    gridsize=25


    counts_collated=read_counts(input_folder=input_folder)
    counts_collated, dict_ratios, CLIC_complexes_count=calc_ratios(counts_collated, pair)
    cm=create_custom_colourmap(cmap, n_colors)

    fig, axes = plt.subplots(nrows, ncols, figsize=(5, 35),
                             #gridspec kwarg allows for the last plot to be smaller, i.e. the colour map colorbar!
                              gridspec_kw={'height_ratios': [4, 4, 4, 4, 4, 4, 0.5], })
    fig, axes=plot_hexbin_kde(fig, axes, counts_collated, vmin, vmax, xlimo, ylimo, cm, cmap_axes, gridsize)
    fig.savefig(f'{output_folder}{pair}.svg')
    plt.show()
