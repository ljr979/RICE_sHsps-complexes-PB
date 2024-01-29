"""This script finds the matched molecules, and the step size from py4bleaching, and finds the size of each matched molecule for plotting. 

This is appropriate only for those molecules that are IN COMPLEXES. i.e., colocalised hsps and clients, which you want to have 'matched' prior to any analysis. 
This script matches their molecule names, and then calculates the number of subunits per spot for those matched molecules, according to the median size calculated in py4bleaching analysis. 

Importantly, this matches based on the IMAGEJ output and the folder nesting that these macros outputs.

This is performed before collating any data. i.e. was run on each individual experiment to output 'molecule_size_for_plotting' which is used to plot the hexbin plots.

"""
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
import random 
from matplotlib import colors

def grab_trajectories_paths_only(input_folder, timepoint_folders):
    """finds all the files with trajectories within them, the full path so it can be read in in the next function

    Args:
        input_folder (str): input folder path
        timepoint_folders (list): list of folders which should each be a time point

    Returns:
        list: a list for hsp, and a list for client trajectories paths
    """
    hsp_trajectories_paths=[]
    client_trajectories_pahts=[]
    #loop over each timepoint folder containing images and trajectories outputs
    for tp_folder in timepoint_folders:
        subfolders=[folder for folder in os.listdir(f'{input_folder}{tp_folder}/')]
        #filtering to find 'Trajectories' folder which is output from ImageJ
        if 'Trajectories' not in subfolders:
            new_input=f'{input_folder}{tp_folder}/647/'
        else: 
            new_input=f'{input_folder}{tp_folder}'    
                
        filenames = next(os.walk(f'{new_input}'), (None, None, []))[2]
        image_folder_list=[folder for folder in os.listdir(f'{new_input}') if folder not in filenames]
        image_folder_list=[folder for folder in image_folder_list if 'Trajectories' not in folder]
        #wanting to loop over ONLY the image folders to get the trajectories file from within the image folder, not the filtered Trajectories folder so we had to find it in order to ignore it and take everything else
        image_folders.append(image_folder_list)
        #find colocalised molecule trajectories
        yikes=[[f'{root}/{name}' for name in files if '_colocal_traj.csv' in name]for root, dirs, files in os.walk(f'{new_input}/')]
        yikes=[item for sublist in yikes for item in sublist]
        #Ignore the ones it gathered from the trajectories filter folder
        yikes=[item for item in yikes if 'Trajectories' not in item]
        yikes=[item for item in yikes if item not in filenames]
        #split up into client and hsp trajectories
        clients=[path for path in yikes if client in path]
        hsp_trajs=[path for path in yikes if hsp in path]
        hsp_trajectories_paths.append(hsp_trajs)
        client_trajectories_pahts.append(clients)


    return hsp_trajectories_paths, client_trajectories_pahts

def read_trajectories(hsp_trajs, clients):
    """matches trajectories to one another and saves path information and timepoint information for later on. gives them unique mol names

    Args:
        hsp_trajs (list): list of paths to trajectories files for each type of protein
        clients (list): list of paths to trajectories files for each type of protein

    Returns:
        df: df with all trajectories and new mol names
    """
    
    new=[]
    for trajectory in hsp_trajs:
        #loop over hsp trajectories paths and read them in
        hsp_trajectories=pd.read_csv(trajectory)
        # append 'hsp' to the end of each molecule name to distinguish from client molecules, as the rest of the name will be matching
        hsp_trajectories.columns = [str(col) + '_hsp' for col in hsp_trajectories.columns]
        for client_traj in clients:
            #now loop through the client molecule files, read in the files that match the image you have just read in for the hsps, and add 'client' to the end of their molecule name
            if client_traj.split('/')[-2]==trajectory.split('/')[-2]:
                client_trajectories=pd.read_csv(client_traj)
                client_trajectories.columns = [str(col) + '_client' for col in client_trajectories.columns]
        
        
        for item in hsp_trajectories:
            #now that they have both been read in as dataframes, loop over EACH trajectory. 
            trajectory_hsp=hsp_trajectories[[item]]
            #find the molecule name/number of the trajectory
            hspnumber=item.split('_')[1]
            for item1 in client_trajectories:
                #loop over client trajectories and do the same
                client_trajectory=client_trajectories[[item1]]
                clientnumber=item1.split('_')[1]
                #compare the two, and if they are the same (i.e., they are matched/in complex), join them together, and give them a matched identifier and save the path info and the timepoint it came from.
                if hspnumber == clientnumber:
                    newtraj=trajectory_hsp.join(client_trajectory)
                    newtraj=newtraj.T
                    newtraj['path_info']=trajectory
                    newtraj['identifier']=random.randint(1,100000)
                    if '' in trajectory.split('/'):
                        newtraj['timepoint (min)']=trajectory.split('/')[-5].split('_')[-1]
                    else: 
                         newtraj['timepoint (min)']=trajectory.split('/')[-3].split('_')[-1]
                    newtraj=newtraj.reset_index().rename(columns={'index':'molecule_number'})
                    new.append(newtraj)
    new=pd.concat(new)
    #enumerate so they don't get overwritten with the next file if there are similar names
    new['molecule_number']=[f'{molecule_number}_{x}' for x, molecule_number in enumerate(new['molecule_number'])]
    return new

def matched_counts(new, step_sizes_hsp, step_sizes_client):
    """calculates the number of subunits per molecule, using the median step size gathered from py4bleaching

    Args:
        new (df): dataframe with the molecule counts in them
        step_sizes_hsp (df): step sizes read in here, but calculated from py4bleaching form all molecules in this entire experiment
        step_sizes_client (df): step sizes read in here, but calculated from py4bleaching form all molecules in this entire experiment

    Returns:
        _type_: _description_
    """
    
    molecule_counts = []
    #loop over each trajectory at each timepoint
    for (molecule,timepoint), df in new.groupby(['molecule_number', 'timepoint (min)']):
        #do first for hsps
        if 'hsp' in molecule:
            #find the maximum value (mean of first 10 values)
            hsp_max_fluorescence_value = df[timepoint_columns].iloc[:,0:10].values.mean()
            # Calculate average number of molecules by mean fluorescence / step size
            molecule_count_hsp=hsp_max_fluorescence_value/int(step_sizes_hsp['step_size'].values[step_sizes_hsp['step_type']=='last_step'])
            #get the unique identifier assigned previously to this molecule and its matching client
            ID=df['identifier']
            #make a dataframe with the important information
            molecule_counts.append(pd.DataFrame([molecule, hsp_max_fluorescence_value, molecule_count_hsp, int(timepoint), df['path_info'].values, ID[0]]).T)
            
        #do the same for clients
        if 'client' in molecule:   
            client_max_fluorescence_value = df[timepoint_columns].iloc[:,0:10].values.mean()
            molecule_count_client=client_max_fluorescence_value/int(step_sizes_client['step_size'].values[step_sizes_client['step_type']=='last_step'])
            ID=df['identifier'].values.tolist()
            ID=ID[0]
            # Calculate average number of molecules by mean fluorescence / step size
            molecule_counts.append(pd.DataFrame([molecule, client_max_fluorescence_value, molecule_count_client, int(timepoint), df['path_info'].values, ID]).T)

    #combine all this information
    molecule_counts = pd.concat(molecule_counts)
    molecule_counts.columns = ['molecule_number', 'max_fluorescence', 'molecule_count', 'timepoint (min)', 'path_info', 'UNIQUE_ID']
    return molecule_counts


if __name__ == "__main__":    
    #location of trajectories (imageJ output)
    input_folder = '/'
    #output location
    output_folder = 'data/Figures/Figure_2/D-matched_hexbins/'
    #location that the py4bleaching analysis output to (example)
    python_input = 'data/example_python_output/two-colour/1_trajectory_analysis/client/'
    if not os.path.exists(output_folder):
            os.makedirs(output_folder)

    #these can be changed according to the name assigned to the trajectories
    client='client'
    hsp='Hsp'
    #change this to your experiment identifier
    exp_num = 'experiment_example'

    #import the step sizes data that was gathered during trajectory analysis script (for HSP first)
    step_sizes_hsp=pd.read_csv(f'{python_input}hsp/fitting_changepoints/median_steps.csv')
    step_sizes_hsp=step_sizes_hsp.drop([col for col in step_sizes_hsp.columns.tolist() if 'Unnamed: 0' in col],axis=1)

    #now gather the step sizes from the client analysis
    step_sizes_client=pd.read_csv(f'{python_input}client/fitting_changepoints/median_steps.csv')
    step_sizes_client=step_sizes_client.drop([col for col in step_sizes_client.columns.tolist() if 'Unnamed: 0' in col],axis=1)

    #find the timepoint folders that the images and trajectories are nested underneath
    timepoint_folders=[folder for folder in os.listdir(f'{input_folder}')if 'beads' not in folder]

    image_folders=[]

    #find the trajectories
    hsp_trajs, clients = grab_trajectories_paths_only(input_folder=input_folder, timepoint_folders=timepoint_folders)


    hsp_trajs=[item for sublist in hsp_trajs for item in sublist]
    clients=[item for sublist in clients for item in sublist]

    #make a new dataframe which has given the clients and hsp trajectories matched names
    new=read_trajectories(hsp_trajs,clients)

    timepoint_columns = [col for col in new.columns.tolist() if col not in ['molecule_number','path_info','timepoint (min)', 'identifier']]

    #turn timepoints into numerical
    timepoint_map={'t0':0, 't15':15, 't30':30, 't60':60, 't240':240, 't420':420}
    new['timepoint (min)']=new['timepoint (min)'].map(timepoint_map)

    #calculate molecule size within complexes (matched)
    molecule_counts=matched_counts(new, step_sizes_hsp, step_sizes_client)


    #pivot the molecule counts data frame so that everything is grouped correctly for plotting (and in longform for seaborn)
    molecule_counts['molecule_type']= molecule_counts['molecule_number'].str.split('_').str[-2]
    #save this sorted dataframe before pivoting
    molecule_counts.to_csv(f'{output_folder}molecule_counts_sorted.csv')

    molecule_pivot=pd.pivot(molecule_counts,index=['timepoint (min)', 'UNIQUE_ID'], columns='molecule_type', values='molecule_count').reset_index()
    #make sure the mol number is a float
    molecule_pivot['client']=molecule_pivot['client'].astype(float)
    molecule_pivot['hsp']=molecule_pivot['hsp'].astype(float)
    #when the counts were turned into integers, the smaller step ones turned into zeros, so need to make them =1  in order to plot this
    molecule_pivot['hsp']=molecule_pivot['hsp'].mask(molecule_pivot['hsp']==0).fillna(1)
    #save for plotting!
    molecule_pivot.to_csv(f'{output_folder}{exp_num}_molecule_counts_for_plotting.csv')

