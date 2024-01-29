"""Finds the matched molecule sizes, and plots them in a kdeplot

"""

import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
from loguru import logger

def concat(molecule_count_list):
    """concatinate molecule counts dataframes

    Args:
        molecule_count_list (list): list of files of molecule counts

    Returns:
        df: collated molecule counts
    """
    alls=[]
    for x in molecule_count_list:
        df=pd.read_csv(f'{x}')
        alls.append(df)

    alls=pd.concat(alls)
    return alls

def plot_kde(counts_collated, output_folder, pair, xlimo, ylimo, ncols, cmap):
    """plots the subunit counts (matched) as a kde (heatmap)

    Args:
        counts_collated (df): collated molecule counts dataframe
        output_folder (str): where to save to
        pair (str): which sHsp + client pair you're plotting
        xlimo (int): x limit for your graph
        ylimo (int): y limit
        ncols (int): how many timepoints you ahve (these are plotted as separate heatmaps)
        cmap (str): colour to plot 
    """
    fig, axes = plt.subplots(ncols=ncols, sharey=True, figsize=(18, 6))
    for x, (timepoint, df) in enumerate(counts_collated.groupby('timepoint (min)')):
        sns.kdeplot(ax=axes[x], x=df["client"],
                    y=df["hsp"], shade=True, cmap=cmap)
        axes[x].set_xlabel("", fontsize=2)
        axes[x].set_ylabel("sHsp", fontsize=24)
        stepsize = 10
        end = xlimo+10
        axes[x].xaxis.set_ticks(np.arange(0, end, stepsize))
        axes[x].set(xlim=(0, xlimo), ylim=(0, ylimo))

        axes[x].set_title(f'{timepoint} (min)', fontsize=20)
        axes[x].tick_params(axis='both', labelsize=20)

    fig.suptitle(
        f'Stoichiometries of {pair} in complex over time', y=1, fontsize=26)
    fig.supxlabel('Client', fontsize=24)
    fig.tight_layout()
    # plt.xlabel('xlabel', fontsize=18)

    plt.savefig(f'{output_folder}kde_plots_timepoints_{pair}.png')

    plt.show()

def truncate_colormap(cmap, minval=0.1, maxval=0.7, n=50):
    """truncates the colourmap so that the bins that are less populated are not white (shortens the range of the colourmap)

    Args:
        cmap (str): colourmap name
        minval (float, optional): bottom cutoff parameter. Defaults to 0.1.
        maxval (float, optional): top cutoff parameter. Defaults to 0.7.
        n (int, optional): bins. Defaults to 50.

    Returns:
        _type_: new, truncated colormap
    """
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

if __name__ == "__main__":    
    input_folder='python_output/1_trajectory_analysis/'
    output_folder='src/two-channels/3_stoichiometries/'

    #find all molecule_counts files in your input folder (this should be where the py4bleaching output was)
    molecule_count_list = [[f'{root}/{filename}' for filename in files if ' molecule_counts.csv' in filename] for root, dirs, files in os.walk(f'{input_folder}')]

    molecule_count_list=[item for sublist in molecule_count_list for item in sublist]
    molecule_counts=concat(molecule_count_list)
    molecule_counts.to_csv(f'{output_folder}combined_mol_counts.csv')


    molecule_pivot=pd.pivot(molecule_counts,index=['treatment', 'Unique_ID'], columns='protein', values='last_step_mol_count').reset_index()
    molecule_pivot['client']=molecule_pivot['client'].astype(float)
    molecule_pivot['hsp']=molecule_pivot['hsp'].astype(float)
    molecule_pivot.to_csv(f'{output_folder}molecule_counts_plots_collated.csv')

    #-------------------
    #PLOTTING
    #E.G fluc and hsp27, but wahtever you want this to be saved as
    plot_kde(counts_collated=molecule_pivot, output_folder=output_folder, pair='FLUC and hsp27', xlimo=60, ylimo=50, ncols=6, cmap=truncate_colormap('Purples'))
