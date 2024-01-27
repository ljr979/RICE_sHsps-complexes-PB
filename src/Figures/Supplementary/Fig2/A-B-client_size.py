"""this script plots the data in Supplementary Figure 2: the difference in client size between the start and end of the incubation + - aB-c (this is panel A-B) AS VIOLINPLOTS.

"""
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger


def filter_for_client_aBc(input_folder, files, grouping_dict):
        
    CLIC=[]
    
    FLUC=[]
    
    for f in files: 
        startfin=pd.read_csv(f'{input_folder}{f}')
        startfin.drop([col for col in startfin.columns.tolist() if 'Unnamed: 0' in col],axis=1, inplace=True)

        startfin=startfin.rename(columns={'pair':'Pair'})
        keeps=[item for item in startfin['Pair'] if 'hsp27' not in item]
        keeps= list(set(keeps))
        keeps = keeps + ['None']
        startfin_abs=startfin[startfin['Pair'].isin(keeps)]
        #map on groups so that we keep coloc and non coloc together, and just label them + sHsp, and control is now -sHsp
        startfin_abs['incubated']=startfin_abs['Colocalisation'].map(grouping_dict)

        #split up all the clients that have been incubated with aB-c
        df_clic=startfin_abs[startfin_abs['Protein']=='CLIC']
        
        df_fluc=startfin_abs[startfin_abs['Protein']=='FLUC']
        

        CLIC.append(df_clic)
        
        FLUC.append(df_fluc)

    CLIC=pd.concat(CLIC)
    
    FLUC=pd.concat(FLUC)
    return CLIC, FLUC

def filter_outlier(df, threshold_dict):
    protein=df['Protein'].tolist()[0]
    threshold=threshold_dict[protein]
    df = df[df['last_step_mol_count'] < threshold]
    return df

def plotting(df, protein, palette, output_folder, chap):
    dfmelt=pd.melt(df, id_vars=['Timepoint','Protein', 'Colocalisation', 'Molecule_number', 'Pair', 'incubated', 'start_end', 'last_step_mol_count'], value_vars=['log_count'])

    fig, ax = plt.subplots()
    ax = sns.violinplot(x="incubated",
        y="value", 
        hue="start_end",
        hue_order=['start', 'end'],
        data=dfmelt,  
        order=['-sHsp','+sHsp'],
        scale='width', 
        palette=palette,
        saturation=0.5)
    ax.set_ylabel('# of subunits')
    ax.set_xlabel(f'{protein} in the absence or presence of {chap}')
    ax.set_ylim(0,3)
    plt.legend(title=f'Incubation time (h)', loc='upper right')
    plt.title(f'{protein} + aBc')

    #plt.savefig(f'{output_folder}_{protein}_log_aBc_stoichiometries.svg')
    plt.show()


if __name__ == "__main__":
    input_folder= 'data/Figures/violinplots-supp-figs2-3-4-coloc-non-coloc-combined/start_finish_filtered/'
    output_folder='data/Figures/violinplots-supp-figs2-3-4-coloc-non-coloc-combined/'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    files=[item for item in os.listdir(input_folder)]
    grouping_dict={'Non-coloc': '+sHsp', 'Coloc':'+sHsp', 'Control':'-sHsp'}
    CLIC, FLUC = filter_for_client_aBc(input_folder, files, grouping_dict)
    # CLIC.to_csv(f'{output_folder}CLIC.csv')
    # FLUC.to_csv(f'{output_folder}FLUC.csv')
    chap='aB-c'
    clients=[CLIC, FLUC]
    colour_dict={'CLIC':'Reds', 'FLUC':'Purples'}

    for b in clients:
        prot_name=b['Protein'].unique()[0]
        palette=colour_dict[prot_name]
        b['log_count']=np.log10(b['last_step_mol_count'])
        plotting(b, prot_name, palette, output_folder)
