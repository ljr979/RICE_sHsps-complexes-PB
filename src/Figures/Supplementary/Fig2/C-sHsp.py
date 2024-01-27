"""plots violinplots of aB-c size in the presence and absence of each of the client proteins, at the start and end of the incubation. Data for panel C in Supp figure 2. 
The difference with this script is that whilst the other ones for aBc and hsp27 compare the client size with and without sHsps, this one will compare the different sHsp sizes with and without clients present. 
"""
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

def filter_for_sHsp_aBc(input_folder, files, grouping_dict):
        
    aBc = []
    
    for f in files: 
        startfin=pd.read_csv(f'{input_folder}{f}')
        startfin.drop([col for col in startfin.columns.tolist() if 'Unnamed: 0' in col],axis=1, inplace=True)

        startfin=startfin.rename(columns={'pair':'Pair'})
        keeps=[item for item in startfin['Pair'] if 'hsp27' not in item]
        keeps= list(set(keeps))
        keeps = keeps + ['None']
        startfin_abs=startfin[startfin['Pair'].isin(keeps)]


        #split up all the clients that have been incubated with aB-c
        df_abc=startfin_abs[startfin_abs['Protein']=='aB-c']
        df_abc2=startfin_abs[startfin_abs['Protein']=='aBc']
        #df_fluc=startfin_abs[startfin_abs['Protein']=='FLUC']
        
        aBc.append(df_abc2)
        aBc.append(df_abc)
        
        

    aBc=pd.concat(aBc)
    #map on groups so that we keep coloc and non coloc together, and just label them + sHsp, and control is now -sHsp
    
    fix_wrong_abc={'CLIC_aB-c':'CLIC_aBc', 'FLUC_aB-c':'FLUC_aBc', 'None':'None','FLUC_aBc':'FLUC_aBc', 'CLIC_aBc':'CLIC_aBc'}
    aBc['Pair']=aBc['Pair'].map(fix_wrong_abc)
    aBc['incubated']=aBc['Pair'].map(grouping_dict)
    #FLUC=pd.concat(FLUC)
    return aBc

def filter_outlier(df, threshold_dict):
    protein=df['Protein'].tolist()[0]
    threshold=threshold_dict[protein]
    df = df[df['last_step_mol_count'] < threshold]
    return df

def plotting2(df,chap):
    dfmelt=pd.melt(df, id_vars=['Timepoint','Protein', 'Colocalisation', 'Molecule_number', 'Pair', 'incubated', 'start_end', '+','last_step_mol_count'], value_vars=['log_count'])

    fig, ax = plt.subplots()
    ax = sns.violinplot(x="+",
        y="value", 
        hue="start_end",
        hue_order=['start', 'end'],
        data=dfmelt,  
        order=['-client','+ CLIC', '+ FLUC'],
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

    chap='aBc'

    files=[item for item in os.listdir(input_folder) ]
    grouping_dict={'FLUC_aBc': '+client', 'CLIC_aBc':'+client', 'None':'-client'}
    aBc = filter_for_sHsp_aBc(input_folder, files, grouping_dict)
    #aBc.to_csv(f'{output_folder}aBc.csv')

    pairs=['CLIC_aBc', 'FLUC_aBc']

    new=[]
    for pair in pairs:
        pair
        p=['None', pair]
        client=pair.split('_')[0]
        new_df=aBc[aBc['Pair'].isin(p)]
        col_dict={'-client':'-client', '+client':f'+ {client}'}
        new_df['+']=new_df['incubated'].map(col_dict)
        new_df['log_count']=np.log10(new_df['last_step_mol_count'])
        new.append(new_df)
    all_aBc=pd.concat(new)

    all_aBc = all_aBc.apply(lambda x: x.apply(lambda y: np.nan if y < 0 else y) if np.issubdtype(x.dtype, np.number) else x)
    all_aBc=all_aBc.dropna(axis=0)

    df=all_aBc

    plotting2(df, chap)