"""This script finds the summary data (mean, sem) and then plots, the colocalisation of each client with aB-c (Panel B, figure 2)

"""
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
from statistics import mean as mean
from scipy.stats import sem

def read_wrangle(input_folder):
    """this function finds and organises all colocalisation data for plotting
    Args:
        input_folder (str): folder containing the colocalisation data for a specific pair of proteins

    Returns:
        df: dataframe with the colocalisation data
    """
    #find all files for this pair of proteins
    coloc_files = [filename for filename in os.listdir(f'{input_folder}') if 'percent_colocalisation_all' in filename]
    coloc_all = []
    #read in and concatinate them
    for filename in coloc_files:
        count = pd.read_csv(f'{input_folder}{filename}')
        coloc_all.append(count)
    coloc_all = pd.concat(coloc_all)
    coloc_all.drop([col for col in coloc_all.columns.tolist() if 'Unnamed: 0' in col],axis=1, inplace=True)
    #rename the timepoints to integers so that they will plot as a lineplot, not categorical
    timepoint_dict = {'zero':0, '15min':0.25, '20min':0.25, '30min':0.5, '40min':0.5,'60min':1,'4h':4, '7h':7, '240min': 4,'420min':7}
    coloc_all['timepoint'] = coloc_all['timepoint'].map(timepoint_dict)
    
    test_client = coloc_all[coloc_all['protein']=='client']
    test = pd.melt(test_client, id_vars=['timepoint','protein'], value_vars='percent_colocalisation')
    return test

def summaries(test): 
    """Find the average between the replicates, and SEM

    Args:
        test (df): dataframe with the colocalisation data

    Returns:
        df: dataframe with the summaries in it
    """
    
    agg_func_math = {
    'value': [ 'mean',  'sem']
    }
    summary_df = test.groupby(['timepoint'], as_index=False).agg(agg_func_math).round(2).reset_index()
    summary_df.columns = ['_'.join(col).rstrip('_') for col in summary_df.columns.values]
    return summary_df

if __name__ == "__main__":
    #define the path to where each of the colocalisation data for each pair lives
    FLUC_aBc = 'data/Figures/Figure_2/B-colocalisation/FLUC_aB-c/'

    CLIC_aBc = 'data/Figures/Figure_2/B-colocalisation/CLIC_aB-c/'
    #where will you save these plots
    output_folder = 'data/Figures/Figure_2/B-colocalisation/'

    #create output folder if it doesn't exist yet
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    #CLIC
    #define the pair you want to work on first (this path to the python results should be defined above)
    test_CLIC = read_wrangle(input_folder=CLIC_aBc)
    #the df to work on here should be the variable you just created in the previous line
    summary_df_CLIC = summaries(test=test_CLIC)
    #removing timepoint 7h because FLUC does not have a 7 h timepoint
    summary_df_CLIC_no7 = summary_df_CLIC[summary_df_CLIC['timepoint']!=7.00]


    #FLUC
    test_FLUC = read_wrangle(input_folder=FLUC_aBc)
    summary_df_FLUC = summaries(test=test_FLUC)
    summary_df_FLUC.to_csv(f'{output_folder}FLUC_aB-c_coloc.csv')

    #PLOTTING
    #colour to use for each protein combination
    colors = ['orange', 'darkviolet']
    #plotting
    fig, ax = plt.subplots()
    #plot CLIC and aB-c points and error
    ax = plt.scatter(x=summary_df_CLIC_no7['timepoint'], y=summary_df_CLIC_no7['value_mean'], c=colors[0], alpha=0.3)
    #plot the error bar at each point
    ax = plt.errorbar(x=summary_df_CLIC_no7['timepoint'], y=summary_df_CLIC_no7['value_mean'], yerr=summary_df_CLIC_no7['value_sem'], c=colors[0], alpha=0.3)

    #now plot FLUC and aB-c points and error
    ax1 = plt.scatter(x=summary_df_FLUC['timepoint'], y=summary_df_FLUC['value_mean'], c=colors[1], alpha=0.5)
    ax1 = plt.errorbar(x=summary_df_FLUC['timepoint'], y=summary_df_FLUC['value_mean'], yerr=summary_df_FLUC['value_sem'], c=colors[1], alpha=0.5)

    #define legend handles
    plt.legend((ax, ax1,),
            ('CLIC', 'FLUC'),
            scatterpoints=1,
            loc='upper right',
            ncol=1,
            fontsize=12)

    #this plots the line that joins the points together, for each protein pair
    ax = plt.plot(summary_df_CLIC_no7['timepoint'], summary_df_CLIC_no7['value_mean'],c='#EA9648')
    ax = plt.plot(summary_df_FLUC['timepoint'], summary_df_FLUC['value_mean'], c='#9888AA')
    ax = plt.tick_params(axis='both', labelsize=16)
    ax = plt.gca()
    ax.set_ylim([0,100])
    plt.ylabel('Percent client colocalisation (%)', fontsize=16)
    plt.title(f'Percentage of client colocalised with aB-c',fontsize=16)
    plt.xlabel('Time (h)')
    ax = ax.get_figure()
    #save!
    ax.savefig(f'{output_folder}mean_scatter_percent_colocal.svg', bbox_inches='tight')