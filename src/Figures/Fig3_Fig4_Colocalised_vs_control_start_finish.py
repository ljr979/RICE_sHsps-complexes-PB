
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

#this script collates all of the non-colocalised and colocalised and control molecule_counts data from py4bleaching, so that we can then plot them together at the start and end of the incubation, without distinguishing between non-colocalised and co-localised status.
def read_counts_controls(stoich_files_controls):
    """Reads in and organises stoichiometry files for controls (i.e., size of clients and sHsps in the absence of each other). The purpose is to make sure the dataframe from these molecules matches the dataframes from colocalised and non-colocalised molecules, so we can make them into a combined dataframe to plot later.

    Args:
        stoich_files_controls (list): list of all the files which have molecule counts in them

    Returns:
        df: dataframe with controls subunit counts
    """
    #filter for columns you wish to keep
    cols_to_keep=[
        'molecule_number',
        'last_step_mol_count', 
        'treatment', 
        'colocalisation',
        'protein', 

         ]
    #read in files and filter for only control conditions
    alls_controls=[]
    for filename in stoich_files_controls:
        count=pd.read_csv(f'{filename}')
        count.drop([col for col in count.columns.tolist() if 'Unnamed: 0' in col],axis=1, inplace=True)
        count=count[count['colocalisation']=='Control']
        count.drop([col for col in count.columns.tolist() if col not in cols_to_keep],axis=1, inplace=True)

        alls_controls.append(count)
    alls_controls=pd.concat(alls_controls)

    #add a column that matches the 'time' that the measurement was taken for easy filtering later
    test=[]
    for row, df in alls_controls.groupby('treatment'):
        row
        df
        timepoint=row.split('-')[0]
        df['timepoint']=timepoint
        test.append(df)
    test=pd.concat(test)

    #turn categorical timepoints into numeric
    tp_dict={'zero': 0,'15min':15, '30min':30, '20min':15, '40min':30, '60min': 60, '180min':180, '240min':180, '420min':420, '4h':180, '7h':420}
    test['timepoint']=test['timepoint'].map(tp_dict)
    test.drop([col for col in test.columns.tolist() if 'treatment' in col],axis=1, inplace=True)

    #filter for the first and last times only
    last_tp=test['timepoint'].unique().max()
    first_tp=test['timepoint'].unique().min()
    times_filter=[first_tp, last_tp]
    #take this subset of df and make a new df
    controls_start_fin=test[test['timepoint'].isin(times_filter)]
    all_start_fin_final=[]
    #add a new column which says which time this was (as start or end, rather than specific time in h)
    for time, df in controls_start_fin.groupby('timepoint'):
        if time == first_tp:
            se = 'start'
        if time == last_tp:
            se ='end'
        df['start_end']=se
        all_start_fin_final.append(df)

    all_start_fin_final=pd.concat(all_start_fin_final)
    #add column to say this was not a 'pair' of protein just a single protein
    all_start_fin_final['Pair']='None'
    #assign a unique molecule name
    test['UNIQUE_ID']= test['molecule_number']=[f'{molecule_number}_{x}' for x, molecule_number in enumerate(test['molecule_number'])]
    #make column names match those from colocalised datasets
    all_start_fin_final=all_start_fin_final.rename(columns={'molecule_number':'Molecule_number', 'timepoint': 'Timepoint', 'protein':'Protein', 'colocalisation':'Colocalisation', 'pair':'Pair'})
    return all_start_fin_final

def read_counts_noncoloc(stoich_files_noncoloc):

    cols_to_keep=[
        'molecule_number',
        'last_step_mol_count', 
        'timepoint', 
        'colocalisation',
        'protein', 

        ]
    noncolocal_counts=[]
    for filename in stoich_files_noncoloc:
        count=pd.read_csv(f'{filename}')
        filename=filename.replace('\\', '/')
        count.drop([col for col in count.columns.tolist() if 'Unnamed: 0' in col],axis=1, inplace=True)
        #leaving this line in incase any pesky colocals infiltrated
        count=count[count['colocalisation']=='Non-coloc']
        #keep only the columns we want to match later
        count.drop([col for col in count.columns.tolist() if col not in cols_to_keep],axis=1, inplace=True)
        #different experiments kind of had different timepoints and ways of saying them i.e. some are numbers some are strings etc so just standardising here
        if type(count['timepoint'].values.tolist()[0])==str:
    
            tp_dict={'zero': 0,'15min':15, '30min':30, '20min':15, '40min':30, '60min': 60, '180min':240, '240min':240, '420min':420, '4h':180, '7h':420}
            
            count['timepoint']=count['timepoint'].map(tp_dict)
        #need to say which PAIR the non coloc came from, because although not colocalised, they are still contributing to the average size of the protein when in the PRESENCE of the other protein in the pair.
        pair=filename.split('/')[-2]
        
        count['pair']=pair

        #save thiso ne
        count.to_csv(f'{output_folder}noncoloc_{pair}_wrangled.csv')
        #find the last timepoint, becasue some experiments couldn't see any molecules at 7h so the last timepoint is different.
        last_tp=count['timepoint'].unique().max()
        first_tp=count['timepoint'].unique().min()
        #filter the dataframe so that it only gives us mol size at these times
        times_filter=[first_tp, last_tp]
        noncoloc_start_fin=count[count['timepoint'].isin(times_filter)]
        #this little loop accounts for the fact that some end points are different! so we just have an extra column that I can group by when plotting, so I know which spot is tart and which is end rather than the actual TIME which will confuse things later
        all_start_fin_final=[]
        for time, df in noncoloc_start_fin.groupby('timepoint'):
            if time == first_tp:
                se = 'start'
            if time == last_tp:
                se ='end'
            df['start_end']=se
            all_start_fin_final.append(df)

        all_start_fin_final=pd.concat(all_start_fin_final)

        noncolocal_counts.append(all_start_fin_final)

    noncolocal_counts=pd.concat(noncolocal_counts)

    noncolocal_counts=noncolocal_counts.rename(columns={'molecule_number':'Molecule_number', 'timepoint': 'Timepoint', 'protein':'Protein', 'colocalisation':'Colocalisation', 'pair':'Pair'})
        
    return noncolocal_counts



    #now need to read in all of the dataframes where I used the matched 

def read_colocal_counts_together(stoich_files_coloc):
    alls_coloc_filtered_start_fin=[]
    for filename in stoich_files_coloc:
        #read in file
        count=pd.read_csv(f'{filename}')
        count.drop([col for col in count.columns.tolist() if 'Unnamed: 0' in col],axis=1, inplace=True)
        #for these colocalised ones, so we can match them to the other dataframes, we just need to melt them and then make sure all the column names match! that is, the 'cliient or hsp' columns need to instead say which protein it is looking at, so we will match that here, for each protein combo
        test_melt=pd.melt(count, 
        id_vars = ['UNIQUE_ID', 'timepoint (min)'],
        value_vars = ['client', 'hsp'], 
        var_name = 'Protein', 
        value_name = 'last_step_mol_count')
        pair = filename.split('/')[-2].split('\\')[-1]
        client=pair.split('_')[0]
        hsp=pair.split('_')[1]
        pair_dict={'client':client, 'hsp':hsp}
        test_melt['Protein']=test_melt['Protein'].map(pair_dict)
        test_melt=test_melt.rename(columns={'UNIQUE_ID':'Molecule_number', 'timepoint (min)': 'Timepoint'})
        test_melt['Colocalisation']='Coloc'
        Exp_num=filename.split('/')[-1].split('_')[0]
        #now all of the columns match the dataframes that have come straight from the py4bleaching analysis!
        #save this total rejigged dataframe for each pair, before filtering them for the first and last timepoint
        test_melt.to_csv(f'{output_folder}{Exp_num}_{pair}_colocal_all_melted.csv')
        last_tp=test_melt['Timepoint'].unique().max()
        first_tp=test_melt['Timepoint'].unique().min()
        times_filter=[first_tp, last_tp]
        
        all_start_fin=test_melt[test_melt['Timepoint'].isin(times_filter)]
        all_start_fin_final=[]
        for time, df in all_start_fin.groupby('Timepoint'):
            if time == first_tp:
                se = 'start'
            if time == last_tp:
                se ='end'
            df['start_end']=se
            all_start_fin_final.append(df)
        all_start_fin_final=pd.concat(all_start_fin_final)
        all_start_fin_final['pair']=pair

        alls_coloc_filtered_start_fin.append(all_start_fin_final)
    alls_coloc_filtered_start_fin=pd.concat(alls_coloc_filtered_start_fin)
    alls_coloc_filtered_start_fin.to_csv(f'{output_folder}colocal_start_fin.csv')
    return alls_coloc_filtered_start_fin

if __name__ == "__main__":    
    input_folder= 'data/'
    output_folder='python_results/coloc_vs_control_start_finsih/coloc_vs_control/'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    stoich_files_controls =[[f'{root}/{filename}' for filename in files if 'molecule_counts.csv' in filename] for root, dirs, files in os.walk(f'{input_folder}')]
    stoich_files_controls=[item for sublist in stoich_files_controls for item in sublist]


    #now filter all the controls so that I just have start and end point!
    controls_only_start_fin=read_counts_controls(stoich_files_controls)
    #save them fore filtering and plotting later
    controls_only_start_fin.to_csv(f'{output_folder}controls_start_fin.csv')


    #----------------------------------------
    #now to find all the ones that are not CONTROL, but are non-coloc and coloc (incubated with each other)
    stoich_files_coloc =[[f'{root}/{filename}' for filename in files if 'molecule_counts_for_plotting.csv' in filename] for root, dirs, files in os.walk(f'{input_folder}')]
    stoich_files_coloc=[item for sublist in stoich_files_coloc for item in sublist]

    colocalised_counts=read_colocal_counts_together(stoich_files_coloc)

    #last step!! need to read in the non-colocalised (but incubated together) counts, and make sure the data structure and column names etc. match, so we can smoosh these all togheter and they will match for plotting. 

    stoich_files_noncoloc =[[f'{root}/{filename}' for filename in files if 'non-coloc_counts.csv' in filename] for root, dirs, files in os.walk(f'{input_folder}')]
    stoich_files_noncoloc=[item for sublist in stoich_files_noncoloc for item in sublist]

    stoich_files_noncoloc=[filename for filename in stoich_files_noncoloc if not 'collated' in filename]

    noncolocal_counts=read_counts_noncoloc(stoich_files_noncoloc)
    noncolocal_counts.to_csv(f'{output_folder}noncoloc_start_fin.csv')