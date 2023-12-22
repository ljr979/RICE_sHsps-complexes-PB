
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
#this script just organises the data for plotting in the next plots. gathers the data which was cleaned and organised in preparation for figures2 and 3, and then combines non colocalised and colocalised files but with separate names so we can split them, and makes sure ALL timepoints are here (in figs 2 and 3, we only had the start and end of the experiment)

input_folder= 'data/Figures/violinplots-supp-figs2-3-coloc-non-coloc-combined/'
output_folder='python_results/Figure_4/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)


stoich_files_noncoloc =[[f'{root}/{filename}' for filename in files if 'wrangled.csv' in filename] for root, dirs, files in os.walk(f'{input_folder}')]
stoich_files_noncoloc=[item for sublist in stoich_files_noncoloc for item in sublist]
stoich_files_coloc =[[f'{root}/{filename}' for filename in files if 'colocal_all_melted.csv' in filename] for root, dirs, files in os.walk(f'{input_folder}')]
stoich_files_coloc=[item for sublist in stoich_files_coloc for item in sublist ]

#this line adds the two lists together so that all the filenames that I outputted from 'fig3_fig4_ colocalisedvscontrol start finish' script are in one big list. THese files are all in the same format, so I should be able to read them in and then fix and filter them so that they can be plotted the way I want them! 
collated_stoich_files=stoich_files_noncoloc+stoich_files_coloc


#Experiment_number='Experiment_65_4'
def read_counts_all(collated_stoich_files):

    alls=[]
    for filename in collated_stoich_files:
        count=pd.read_csv(f'{filename}')
        filename=filename.replace('//', '/')
        info = filename.split('/')[-1]
        count.drop([col for col in count.columns.tolist() if 'Unnamed: 0' in col],axis=1, inplace=True)
        #now need to make a dictionary so that all the columns are named the same!!! otherwise they can't be concatinated properly- some seem to have capital letters etc, but need to match completely
        proteins=info.split('_')[1:3]
        pair=proteins[0]+'_'+proteins[1]
        count['pair']=pair
        count=count.rename(columns={'molecule_number':'Molecule_number','timepoint':'Timepoint','protein':'Protein','colocalisation':'Colocalisation', 'pair':'Pair'})


        #count.drop([col for col in count.columns.tolist() if col not in cols_to_keep],axis=1, inplace=True)

        alls.append(count)
    alls=pd.concat(alls)

    tp_dict={0: 0, 15:20, 20:20, 30:40, 40:40, 60: 60, 180:240, 240:240, 420:420}
    
    alls['Timepoint']=alls['Timepoint'].map(tp_dict)
    pairs_dict={
        'CLIC_aBc':'CLIC_aB-c',
        'FLUC_aBc': 'FLUC_aB-c', 
        'CLIC_hsp27':'CLIC_hsp27',
        'FLUC_hsp27':'FLUC_hsp27',
        'Rhod_hsp27':'Rhod_hsp27',
        'CLIC_aB-c':'CLIC_aB-c',
        'FLUC_aB-c':'FLUC_aB-c'}
    alls['Pair']=alls['Pair'].map(pairs_dict)
    return alls

all_proteins_tps=read_counts_all(collated_stoich_files)


#save them fore filtering and plotting later
all_proteins_tps.to_csv(f'{output_folder}allproteins_alltimepoints_coloc_noncoloc.csv')



