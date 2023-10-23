       
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns
from loguru import logger

input_folder='python_output/1_trajectory_analysis/'
output_folder='src/two-channels/3_stoichiometries/'

molecule_count_list = [[f'{root}/{filename}' for filename in files if ' molecule_counts.csv' in filename] for root, dirs, files in os.walk(f'{input_folder}')]

molecule_count_list=[item for sublist in molecule_count_list for item in sublist]

def concat(molecule_count_list):


    alls=[]
    for x in molecule_count_list:
        df=pd.read_csv(f'{x}')
        alls.append(df)

    alls=pd.concat(alls)
    return alls


def plot_kde(counts_collated, output_folder, pair, xlimo, ylimo, ncols, cmap):
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
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

molecule_counts=concat(molecule_count_list)
molecule_counts.to_csv(f'{output_folder}combined_mol_counts.csv')


molecule_pivot=pd.pivot(molecule_counts,index=['treatment', 'Unique_ID'], columns='protein', values='last_step_mol_count').reset_index()
molecule_pivot['client']=molecule_pivot['client'].astype(float)
molecule_pivot['hsp']=molecule_pivot['hsp'].astype(float)
molecule_pivot.to_csv(f'{output_folder}molecule_counts_for_plotting.csv')


#-------------------
#PLOTTING
#E.G fluc and hsp27, but wahtever you want this to be saved as
pair = 'FLUC and hsp27'

#whatever you want also.
cmap = 'Purples'

#can change these to zoom in etc.
ylimo = 50
xlimo = 60
ncols = 6
counts_collated = molecule_pivot
cmap = plt.get_cmap(cmap)
new_cmap = truncate_colormap(cmap)
plot_kde(counts_collated, output_folder, pair, xlimo, ylimo, ncols, cmap)
