"""This plots, for each client incubated with hsp27, the # Of molecules  +/- Hsp27  at  7 h incubation for all combined controls and complexes and non-colocalised molecules (Figure 3, Panel B)

"""
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

def filter_for_client_hsp27(input_folder, files, grouping_dict):
    """filter the Hsp27 dataframes to only look at the CLIENT rather than the hsp. return these dataframes for plotting

    Args:
        input_folder (str): folder where these files live
        files (list): files to read in and combine
        grouping_dict (dict): a dictionary which tells this function how to rename the data for plotting

    Returns:
        df: 3 dfs which have the data read to plot
    """
    CLIC=[]
    FLUC=[]
    RHOD=[]
    for f in files: 
        startfin=pd.read_csv(f'{input_folder}{f}')
        startfin.drop([col for col in startfin.columns.tolist() if 'Unnamed: 0' in col],axis=1, inplace=True)

        startfin=startfin.rename(columns={'pair':'Pair'})
        #don't want to include aBc here just hsp27
        keeps=[item for item in startfin['Pair'] if 'aB-c' not in item]
        keeps= list(set(keeps))
        keeps = keeps + ['None']
        startfin_abs=startfin[startfin['Pair'].isin(keeps)]
        #map on groups so that we keep coloc and non coloc together, and just label them + sHsp, and control is now -sHsp
        startfin_abs['incubated']=startfin_abs['Colocalisation'].map(grouping_dict)

        #split up all the clients that have been incubated with Hsp27
        df_clic=startfin_abs[startfin_abs['Protein']=='CLIC']
        df_rhod=startfin_abs[startfin_abs['Protein']=='Rhodanese']
        df_fluc=startfin_abs[startfin_abs['Protein']=='FLUC']
        

        CLIC.append(df_clic)
        RHOD.append(df_rhod)
        FLUC.append(df_fluc)

    CLIC=pd.concat(CLIC)
    RHOD=pd.concat(RHOD)
    FLUC=pd.concat(FLUC)
    return CLIC, FLUC, RHOD

def plotting(df, protein, palette, output_folder):
    """plot as violinplots the log of subunit count

    Args:
        df (df): dataframe with the combined molecule counts for each client
        protein (str): the client protein which is being plotted
        palette (str): colour to plot
        output_folder (str): where to save
    """
    dfmelt=pd.melt(df, id_vars=['Timepoint','Protein', 'Colocalisation', 'Molecule_number', 'Pair', 'incubated', 'start_end', 'last_step_mol_count'], value_vars=['log_count'])

    fig, ax = plt.subplots()
    ax = sns.violinplot(x="Protein",
        y="value", 
        hue="incubated",
        hue_order=['-sHsp', '+sHsp'],
        data=dfmelt,  
        order=['CLIC','FLUC', 'Rhodanese'],
        scale='width', 
        palette=palette,
        saturation=0.5)
    ax.set_ylabel('# of subunits')
    ax.set_xlabel(f'{protein} in the absence or presence of hsp27')
    ax.set_ylim(0,3)
    plt.legend(title=f'Incubation time (h)', loc='upper right')
    plt.title(f' hsp27')

    plt.savefig(f'{output_folder}_{protein}_log_hsp27_stoichiometries.svg')
    plt.show()

if __name__ == "__main__":
    input_folder='data/Figures/violinplots-supp-figs2-3-coloc-non-coloc-combined/start_finish_filtered/'
    output_folder='data/Figures/Figure_3/B-violinplots/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    files=[item for item in os.listdir(input_folder)]
    grouping_dict={'Non-coloc': '+sHsp', 'Coloc':'+sHsp', 'Control':'-sHsp'}
    #filter and organise, and then save these
    CLIC, FLUC , RHOD=filter_for_client_hsp27(input_folder, files, grouping_dict)
    CLIC.to_csv(f'{output_folder}CLIC.csv')
    FLUC.to_csv(f'{output_folder}FLUC.csv')
    RHOD.to_csv(f'{output_folder}RHOD.csv')
    #these lines filter my previous dfs from being CLIC and fluc alone, with start and finish, to being just the END point, and then concatinates them so that clic and fluc are together, and just the end point.
    c2=CLIC[CLIC['start_end']=='end']
    f2=FLUC[FLUC['start_end']=='end']
    r2=RHOD[RHOD['start_end']=='end']
    combo_v2 = pd.concat([c2, f2, r2])
    #convert to the log of the number of subunits for plotting!
    combo_v2['log_count']=np.log10(combo_v2['last_step_mol_count'])
    #save this data that is going to be plotted
    combo_v2.to_csv(f'{output_folder}violinplots-fig3.csv')
    #plot and save
    plotting(df=combo_v2, protein='CLIC-FLUC-rhod-hsp27', palette='Greens', output_folder=output_folder)
