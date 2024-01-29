"""This script finds the summary data (mean, sem) and then plots, the colocalisation of each client with hsp27 (Panel A, Figure 3)

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
    """ this function finds and organises all colocalisation data for plotting


    Args:
        input_folder (str): input for each pair of proteins

    Returns:
        df: dataframe with the collated colocalisation data
    """

    #find all files for this pair of proteins    
    coloc_files=[filename for filename in os.listdir(f'{input_folder}') if 'percent_colocalisation_all' in filename]
    #read in and concatinate
    coloc_all=[]
    for filename in coloc_files:
        count=pd.read_csv(f'{input_folder}{filename}')
        coloc_all.append(count)
    coloc_all=pd.concat(coloc_all)
    coloc_all.drop([col for col in coloc_all.columns.tolist() if 'Unnamed: 0' in col],axis=1, inplace=True)
    #rename timepoints so that they are numerical not categorical
    timepoint_dict= {'zero':0, '15min':0.25, '20min':0.25, '30min':0.5, '40min':0.5,'60min':1,'4h':4, '7h':7, '240min': 4,'420min':7}
    coloc_all['timepoint']=coloc_all['timepoint'].map(timepoint_dict)

    test_client=coloc_all[coloc_all['protein']=='client']
    test=pd.melt(test_client, id_vars=['timepoint','protein'], value_vars='percent_colocalisation')
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
    summary_df=test.groupby(['timepoint'], as_index=False).agg(agg_func_math).round(2).reset_index()
    summary_df.columns = ['_'.join(col).rstrip('_') for col in summary_df.columns.values]
    return summary_df

if __name__ == "__main__":
    #define the INPUT PATH for each of the protein pairs (where does the colocalisation data live for this pair? send the script there). This should be the output from 'colocalisation analysis' script.
    FLUC_hsp27='data/Figures/Figure_3/A-colocalisation/FLUC_hsp27/'

    Rhod_hsp27='data/Figures/Figure_3/A-colocalisation/Rhod_hsp27/'

    CLIC_hsp27='data/Figures/Figure_3/A-colocalisation/CLIC_hsp27/'

    output_folder='data/Figures/Figure_3/A-colocalisation/hsp27/'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)


    #CLIC
    test_CLIC=read_wrangle(input_folder=CLIC_hsp27)
    summary_df_CLIC=summaries(test=test_CLIC)

    #FLUC
    test_FLUC=read_wrangle(input_folder=FLUC_hsp27)
    summary_df_FLUC=summaries(test=test_FLUC)

    #Rhod
    test_Rhod=read_wrangle(input_folder=Rhod_hsp27)
    summary_df_Rhod=summaries(test=test_Rhod)

    #colour for each line plot
    colors = ['orange', 'darkviolet', 'green']
    #plotting
    #plot scatterplot for CLIC
    fig, ax = plt.subplots()
    ax=plt.scatter(x=summary_df_CLIC['timepoint'], y=summary_df_CLIC['value_mean'], c=colors[0], alpha=0.3)
    #plot error for CLIC
    ax=plt.errorbar(x=summary_df_CLIC['timepoint'], y=summary_df_CLIC['value_mean'], yerr=summary_df_CLIC['value_sem'], c=colors[0], alpha=0.3)

    #now FLUC scatter and error
    ax1=plt.scatter(x=summary_df_FLUC['timepoint'], y=summary_df_FLUC['value_mean'], c=colors[1], alpha=0.5)
    ax1=plt.errorbar(x=summary_df_FLUC['timepoint'], y=summary_df_FLUC['value_mean'], yerr=summary_df_FLUC['value_sem'], c=colors[1], alpha=0.5)

    #now rhodanese scatter and error
    ax2=plt.scatter(x=summary_df_Rhod['timepoint'], y=summary_df_Rhod['value_mean'], c=colors[2], alpha=0.3)
    ax2=plt.errorbar(x=summary_df_Rhod['timepoint'], y=summary_df_Rhod['value_mean'], yerr=summary_df_Rhod['value_sem'], c=colors[2], alpha=0.5)

    #define the legend handles
    plt.legend((ax, ax1, ax2),
            ('CLIC', 'FLUC', 'Rhodanese'),
            scatterpoints=1,
            loc='upper right',
            
            ncol=1,
            fontsize=12)

    #now plot the line that joins the scatter points for each protein
    ax=plt.plot(summary_df_CLIC['timepoint'], summary_df_CLIC['value_mean'],c='#EA9648')
    ax=plt.plot(summary_df_FLUC['timepoint'], summary_df_FLUC['value_mean'], c='#9888AA')
    ax=plt.plot(summary_df_Rhod['timepoint'], summary_df_Rhod['value_mean'], c='#46B030')
    ax=plt.tick_params(axis='both', labelsize=16)
    ax=plt.gca()
    ax.set_ylim([0,100])
    #labels 
    plt.ylabel('Percent client colocalisation (%)', fontsize=16)
    plt.title(f'Percentage of client colocalised with Hsp27',fontsize=16)
    plt.xlabel('Time (h)', fontsize=16)

    ax=ax.get_figure()
    #save!
    ax.savefig(f'{output_folder}Hsp27_mean_colocalisation.svg', bbox_inches='tight')