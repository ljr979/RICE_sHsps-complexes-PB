

import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
#this script needs to filter all the data for coloc vs non-coloc, in the presence of aB-c, and in the case of CLIC, FLUC AND aB-c, then save the data as with figs 3 and 4 for plotting in violinplots.

input_folder= 'D:/Papers/Complexes_paper/python_results/Fig6_noncolocal_vs_colocal/filtering/'
output_folder='python_results/Fig6_noncolocal_vs_colocal/filtering/abc/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

files=[item for item in os.listdir(input_folder) if 'coloc_noncoloc.csv' in item]

def filter_for_sHsp_aBc(input_folder, files):
        
    aBc = []
    
    for f in files: 
        aBc_only=pd.read_csv(f'{input_folder}{f}')
        aBc_only.drop([col for col in aBc_only.columns.tolist() if 'Unnamed: 0' in col],axis=1, inplace=True)

        aBc_only=aBc_only.rename(columns={'pair':'Pair'})
        keeps=[item for item in aBc_only['Pair'] if 'aB-c' in item]
        keeps= list(set(keeps))
        # keeps = keeps + ['None']
        aBc_only_abs=aBc_only[aBc_only['Pair'].isin(keeps)]


        #split up all the clients that have been incubated with hsp27
        #df_aBc=aBc_only_abs[aBc_only_abs['Protein']=='aB-c']
        


        aBc.append(aBc_only_abs)
        
        

    aBc=pd.concat(aBc)
    #map on groups so that we keep coloc and non coloc together, and just label them + sHsp, and control is now -sHsp
    
    
    #FLUC=pd.concat(FLUC)
    return aBc

def filter_outlier(df, threshold_dict):
    protein=df['Protein'].tolist()[0]
    threshold=threshold_dict[protein]
    df = df[df['last_step_mol_count'] < threshold]
    return df

def get_medians_for_plotting(sHsp_only, sHsp):


    sHsp=sHsp_only[sHsp_only['Protein'].isin(sHsp)]
    df1=sHsp

    alls=[]

    for (protein, timepoint, colocalisation), df in df1.groupby(['Pair', 'Timepoint', 'Colocalisation']): 
        protein
        timepoint
        df
        colocalisation
        LIST=[]

        median=df['last_step_mol_count'].median()
        error=df['last_step_mol_count'].sem()
        pair=df['Pair'].unique().tolist()[0]


        LIST=protein, timepoint, colocalisation, pair, median, error
        alls.append(LIST)

    alls=pd.DataFrame(alls, columns=['protein', 'timepoint', 'colocalisation', 'pair', 'median', 'error'])
    return alls

def plot_lines_sHsp(alls, axes_dict, title, output_folder, CLIC_cols, FLUC_cols):
    fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(9,11))
    axes_dict=axes_dict
    for group, df in alls.groupby(['pair']): 

        axo=axes_dict[group]

        group
        
        client=group.split('_')[0]
        sHsp=group.split('_')[1]

        df

        c=df[df['colocalisation']=='Coloc']
        n=df[df['colocalisation']=='Non-coloc']

        if 'CLIC' in group: 
            colours_dict=CLIC_cols  
        if 'FLUC' in group:
            colours_dict=FLUC_cols


        # ax1=sns.lineplot(
        #         data=df[df['colocalisation']=='Coloc'], x="timepoint",
        #         y="median", 
        #         hue='colocalisation',
        #         palette='Greys', ax=axes[axo]
        #         )
        ax1=axes[axo].errorbar(x=c['timepoint'], y=c['median'], yerr=c['error'], c=colours_dict['Coloc'])
        # ax2=sns.lineplot(data=df[df['colocalisation']=='Non-coloc'], x="timepoint",
        #         y="median", 
        #         hue='colocalisation',
        #         palette='Greys', ax=axes[axo]
        #         )
        ax2=axes[axo].errorbar(x=n['timepoint'], y=n['median'], yerr=n['error'], c=colours_dict['Non-coloc'])

        axes[axo].set_ylabel("")
        axes[axo].set_xlabel("")
        axes[axo].xaxis.set_ticks(np.arange(0, 250, 60))
        axes[axo].set(xlim=(0,250), ylim=(0,max(df['median'])+2))
        

        axes[axo].set_title(f'{client} + {sHsp}', y=0.9, va='top')
        axes[axo].tick_params(axis='both')
        axes[axo].legend((ax1, ax2), ('Coloc', 'Non-coloc'),title=f'Colocalisation', loc='upper right', framealpha=0.5 )

            

        fig.suptitle(f'{title}')
        fig.supxlabel('Time (min)')
        fig.supylabel('median # subunits')
        fig.tight_layout()


    plt.savefig(f'{output_folder}{group}_lineplots_subplots_{sHsp}_sHsp.svg')
    plt.savefig(f'{output_folder}{group}_lineplots_subplots_{sHsp}_sHsp.png')

    plt.show()

def get_medians_for_plotting_clients(clients, aBc_only):


    CLIENTS=aBc_only[aBc_only['Protein'].isin(clients)]
    df1=CLIENTS

    alls=[]

    for (protein, timepoint, colocalisation), df in df1.groupby(['Protein', 'Timepoint', 'Colocalisation']): 
        protein
        timepoint
        df
        colocalisation
        LIST=[]

        median=df['last_step_mol_count'].median()
        error=df['last_step_mol_count'].sem()
        pair=df['Pair'].unique().tolist()[0]


        LIST=protein, timepoint, colocalisation, pair, median, error
        alls.append(LIST)

    alls=pd.DataFrame(alls, columns=['protein', 'timepoint', 'colocalisation', 'pair', 'median', 'error'])



    return alls

def plot_lines_clients(alls_clients, axes_dict, title, output_folder, CLIC_cols, FLUC_cols):
    fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(9,11))
    axes_dict=axes_dict
    for group, df in alls_clients.groupby(['pair']): 
        axo=axes_dict[group]
        group
        df
        c=df[df['colocalisation']=='Coloc']
        n=df[df['colocalisation']=='Non-coloc']
        if 'CLIC' in group: 
            colours_dict=CLIC_cols  
        if 'FLUC' in group:
            colours_dict=FLUC_cols
        # ax1=sns.lineplot(
        #         data=df[df['colocalisation']=='Coloc'], x="timepoint",
        #         y="median", 
        #         hue='colocalisation',
        #         palette='Greys', ax=axes[axo]
        #         )
        ax1=axes[axo].errorbar(x=c['timepoint'], y=c['median'], yerr=c['error'], c=colours_dict['Coloc'])
        # ax2=sns.lineplot(data=df[df['colocalisation']=='Non-coloc'], x="timepoint",
        #         y="median", 
        #         hue='colocalisation',
        #         palette='Greys', ax=axes[axo]
        #         )
        ax2=axes[axo].errorbar(x=n['timepoint'], y=n['median'], yerr=n['error'], c=colours_dict['Non-coloc'])


        axes[axo].set_ylabel("")
        axes[axo].set_xlabel("")
        axes[axo].xaxis.set_ticks(np.arange(0, 250, 60))
        axes[axo].set(xlim=(0,250), ylim=(0,max(df['median'])+2))
        

        axes[axo].set_title(f'{group}', y=0.9, va='top')
        axes[axo].tick_params(axis='both')
        axes[axo].legend((ax1, ax2), ('Coloc', 'Non-coloc'),title=f'Colocalisation', loc='upper right', framealpha=0.5 )

            

        fig.suptitle(f'{title}')
        fig.supxlabel('Time (min)')
        fig.supylabel('median # subunits')
        fig.tight_layout()


    plt.savefig(f'{output_folder}{group}_lineplots_subplots_aB-c_clients.svg')
    plt.savefig(f'{output_folder}{group}_lineplots_subplots_aB-c_clients.png')

    plt.show()

#this filters the big dataframe for only those that have been incubated with aB-c
aBc_only = filter_for_sHsp_aBc(input_folder, files)

aBc_only.to_csv(f'{output_folder}aBc_coloc_noncoloc_alltps.csv')

#define the clients present
clients=['FLUC','CLIC']
#determine the median and error for each timepoint, colocalised and non-colocalised, for each protein in these experiments

alls_clients=get_medians_for_plotting_clients(clients, aBc_only)
#made a dictionary to determine which pair you want in which position on your plot- change accordingly! position 0 is top graph, 1 is bottom graph, 2 below that etc. and just change keys to be whatever your pair is defined as
axes_dict={'CLIC_aB-c':0,'FLUC_aB-c':1}
#dictionary to map the client protein colour scheme to the colocalised and non-colocalised conditions within the function
FLUC_cols={'Coloc': '#564787','Non-coloc':'#dbcbd8'}

CLIC_cols={'Coloc': '#fe5d26','Non-coloc':'#f2c078'}  

#title of your graph
title='Median coloc vs non-coloc client (incubated w aB-c) size'
#where you want to save the plot!
output_folder='python_results/Fig6_noncolocal_vs_colocal/filtering/abc/'
#run the plotting function :)
plot_lines_clients(alls_clients, axes_dict, title, output_folder,CLIC_cols,FLUC_cols)



#define the sHsp as it is in the dataframe so you can select for it!
sHsp=['aB-c']
#now we do the same for aB-c 
alls=get_medians_for_plotting(aBc_only, sHsp)
#match this to the clients plot
axes_dict={'CLIC_aB-c':0,'FLUC_aB-c':1}
#title?
title='Median aB-c molecule size, coloc vs. non-coloc'

#plot!
plot_lines_sHsp(alls, axes_dict, title, output_folder,CLIC_cols,FLUC_cols)

