
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger


input_folder= 'data/Figures/violinplots-supp-figs2-3-coloc-non-coloc-combined/start_finish_filtered/'
output_folder='data/Figures/Figure_2/C-violinplots/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)


def filter_for_client_aBc(input_folder, files, grouping_dict):
    """filters counts data to only include aB-c incubated mols

    Args:
        input_folder (str): folder containing filtered/organised counts data
        files (list): files to be read in
        grouping_dict (dict): dictionary to define what to group the incubation categories under

    Returns:
        df: dataframes which are filtered appropriately for plotting, one for each client incubated with aB-c
    """
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


def plotting(df, protein, palette, output_folder):
    #plot as violinplots the log of subunit count
    dfmelt=pd.melt(df, id_vars=['Timepoint','Protein', 'Colocalisation', 'Molecule_number', 'Pair', 'incubated', 'start_end', 'last_step_mol_count'], value_vars=['log_count'])

    fig, ax = plt.subplots()
    ax = sns.violinplot(x="Protein",
        y="value", 
        hue="incubated",
        hue_order=['-sHsp', '+sHsp'],
        data=dfmelt,  
        order=['CLIC','FLUC'],
        scale='width', 
        palette=palette,
        saturation=0.5)
    ax.set_ylabel('# of subunits')
    ax.set_xlabel(f'aB-c')
    ax.set_ylim(0,3)
    plt.legend(title=f'Incubation time (h)', loc='upper right')
    plt.title(f'{protein} + aBc')

    plt.savefig(f'{output_folder}_{protein}_log_aBc_stoichiometries.svg')
    plt.show()

files=[item for item in os.listdir(input_folder)]
grouping_dict={'Non-coloc': '+sHsp', 'Coloc':'+sHsp', 'Control':'-sHsp'}
#filter for clients of interest
CLIC, FLUC = filter_for_client_aBc(input_folder, files, grouping_dict)
CLIC.to_csv(f'{output_folder}CLIC.csv')
FLUC.to_csv(f'{output_folder}FLUC.csv')


#these lines filter my previous dfs from being CLIC and fluc alone, with start and finish, to being just the END point, and then concatinates them so that clic and fluc are together, and just the end point.
c2 = CLIC[CLIC['start_end'] == 'end']
f2 = FLUC[FLUC['start_end'] == 'end']
combo_v2=pd.concat([c2, f2])

#convert count to log10 count for plotting
combo_v2['log_count'] = np.log10(combo_v2['last_step_mol_count'])
#save this file
combo_v2.to_csv(f'{output_folder}violinplots-fig2.csv')
df=combo_v2
palette='Purples'
protein='CLIC-FLUC-aBc'
plotting(df, protein, palette, output_folder)

