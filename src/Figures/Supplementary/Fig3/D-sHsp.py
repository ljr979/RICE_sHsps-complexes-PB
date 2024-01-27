"""script to plot the size of hsp27 molecules in presence and absence of each client, at the start and end of the incubation. This is the data in panel D of Supp figure 3.

"""
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

def filter_for_sHsp_hsp27(input_folder, files, grouping_dict):
    """filters the big dataset with all molecule sizes, so that only those that are hsp27 (not aBc and hsp27) are present

    Args:
        input_folder (str): path to mol sizes filtered for start and end
        files (list): list of files in the input
        grouping_dict (dict): the grouping into + / - client 

    Returns:
        df: dataframe filtered for data only hsp27 related
    """
    hsp27 = []
    
    for f in files: 
        startfin=pd.read_csv(f'{input_folder}{f}')
        startfin.drop([col for col in startfin.columns.tolist() if 'Unnamed: 0' in col],axis=1, inplace=True)

        startfin=startfin.rename(columns={'pair':'Pair'})
        keeps=[item for item in startfin['Pair'] if 'hsp27' in item]
        keeps= list(set(keeps))
        keeps = keeps + ['None']
        startfin_abs=startfin[startfin['Pair'].isin(keeps)]


        #split up all the clients that have been incubated with hsp27
        df_hsp27=startfin_abs[startfin_abs['Protein']=='hsp27']
        


        hsp27.append(df_hsp27)
        
        

    hsp27=pd.concat(hsp27)
    #map on groups so that we keep coloc and non coloc together, and just label them + sHsp, and control is now -sHsp
  
    hsp27['incubated']=hsp27['Pair'].map(grouping_dict)
    #FLUC=pd.concat(FLUC)
    return hsp27

def plotting2(df,chap):
    """plots the molecule sizes of hsp27, as violinplots

    Args:
        df (df): dataframe with the mol sizes
        chap (str): the name of the chaperone, for labelling
    """
    dfmelt=pd.melt(df, id_vars=['Timepoint','Protein', 'Colocalisation', 'Molecule_number', 'Pair', 'incubated', 'start_end', '+','last_step_mol_count'], value_vars=['log_count'])

    fig, ax = plt.subplots()
    ax = sns.violinplot(x="+",
        y="value", 
        hue="start_end",
        hue_order=['start', 'end'],
        data=dfmelt,  
        order=['-client','+ CLIC', '+ FLUC', '+ Rhod'],
        scale='width', 
        palette='Blues',
        saturation=0.5)
    ax.set_ylabel('Log10 (# of subunits)')
    ax.set_xlabel(f'treatment')
    ax.set_ylim(-0.2,max(dfmelt['value']))
    plt.legend(title=f'Incubation time (h)', loc='upper left')
    plt.title(f'{chap} + or - client')

    #plt.savefig(f'{output_folder}hsp27_stoichiometries_plus_minus_client.svg')
    plt.show()

if __name__ == "__main__":

    input_folder= 'data/Figures/violinplots-supp-figs2-3-4/start_finish_filtered/'
    output_folder='data/Figures/violinplots-supp-figs2-3-4/'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    chap='hsp27'

    files=[item for item in os.listdir(input_folder) ]
    grouping_dict={'FLUC_hsp27': '+client', 'CLIC_hsp27':'+client', 'Rhod_hsp27':'+client', 'None':'-client'}
    hsp27 = filter_for_sHsp_hsp27(input_folder, files, grouping_dict)
    pairs=['CLIC_hsp27', 'FLUC_hsp27', 'Rhod_hsp27']

    new=[]
    for pair in pairs:
        pair
        p=['None', pair]
        client=pair.split('_')[0]
        new_df=hsp27[hsp27['Pair'].isin(p)]
        col_dict={'-client':'-client', '+client':f'+ {client}'}
        new_df['+']=new_df['incubated'].map(col_dict)
        new_df['log_count']=np.log10(new_df['last_step_mol_count'])
        new.append(new_df)
    all_hsp27=pd.concat(new)

    all_hsp27 = all_hsp27.apply(lambda x: x.apply(lambda y: np.nan if y < 0 else y) if np.issubdtype(x.dtype, np.number) else x)
    all_hsp27=all_hsp27.dropna(axis=0)

    df=all_hsp27

    plotting2(df, chap)