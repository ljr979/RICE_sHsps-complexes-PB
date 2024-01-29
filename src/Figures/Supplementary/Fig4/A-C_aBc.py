"""This script needs to filter all the data for coloc vs non-coloc, in the presence of aB-c, and in the case of CLIC, FLUC AND aB-c, then save the data as with figs 3 and 4 for plotting in violinplots.
"""
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

def plotting(df1, pair, protein, sHsp, palette_dict):
    """plots the coloc v non-coloc MOL SIZE data as a violinplot , for each protein incubated with aB-c (including aB-c)

    Args:
        df1 (df): dataframe with the mol size information of all molecules incubated with aB-c
        pair (str): the pair incubated to produce these molecule sizes
        protein (str): the protein being plotted
        sHsp (str): the chaperone being incubated with
        palette_dict (dict): the colour to plot the data in
    """
    for group, df in df1.groupby(['Protein']): 
        dfmelt=pd.melt(df, id_vars=['Timepoint','Protein', 'Colocalisation', 'Molecule_number', 'Pair'], value_vars=['last_step_mol_count'])
        group
        pal1=palette_dict[group]
        sns.set_palette(sns.color_palette(pal1))
        palette=pal1
        fig, ax = plt.subplots()
        ax = sns.violinplot(x="Timepoint",
            y="value", 
            hue="Colocalisation",
            #hue_order=['start', 'end'],
            data=dfmelt,  
            #order=['-client','+client'],
            scale='width', 
            palette=palette,
            saturation=0.5)
        ax.set_ylabel('# of subunits')
        ax.set_xlabel(f'Time (min)')
        ax.set_ylim(0,max(dfmelt['value']))
        plt.legend(title=f'{protein} + {sHsp}', loc='upper left')
        plt.title(f'{group}')

        #plt.savefig(f'{output_folder}{pair}_{group}_noncoloc_coloc_stoichiometries.svg')
        plt.show()

if __name__ == "__main__":
    input_folder= 'data/Figures/Figure_4/'

    palette_store={'CLIC_aB-c':['#ff4800', "#ff7900"], 'FLUC_aB-c':['#4b9cd3', '#1b264f'], 'FLUC':'Purples', 'CLIC':'Reds'}
    chap='aB-c'
    df=pd.read_csv(f'{input_folder}{chap}/{chap}_coloc_noncoloc_alltps.csv')
    pairs=[f'FLUC_{chap}',f'CLIC_{chap}']

    for pair in pairs:
        client=pair.split('_')[0]
        df1=df[df['Pair']==pair]
        palette_dict={f'{chap}':palette_store[pair], f'{client}':palette_store[client]}
        plotting(df1, pair, protein=client, sHsp='aB-c', palette_dict=palette_dict)

