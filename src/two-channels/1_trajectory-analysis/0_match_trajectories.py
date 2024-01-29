import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
import random
#this input should be the location of raw trajectory files, for coloc and non-colocalised proteins (both client and sHsp), at each timepoint.(/Trajectories folder content from the imagej processing script). 
#this should be structured as /data/timepoint/coloc&non-coloc/.csv files (these files should be numbered / labelled as the NAME)

def grab_trajectories_paths_only(inputt, client, hsp, tp_folder):
    """Function to grab all trajectory files and filter for those only colocalised.

    Args:
        input_folder (str): input folder where RAW unprocessed trajectories (output from imageJ processing analysis), this should exist for each timepoint (or treatment) in the experiment, and within each of those folders should be a Trajectories folder containing subfolders for client and sHsp, and colocalised and not colocalised proteins trajectories.
        timepoint_folders (list): list of all the 'treatments' folders in this experiment (e.g. timepoints), so we can interate over them to get all the files
        client (string): string defined at start of script which is the name of the FOLDER within colocalisation folder, which contains your client trajectories (this string may change if you change it in the imageJ script, so it is easily changeable in this way)
        hsp (string): see above.
    Returns:
        lists: a list of hsp trajectories, and client trajectories
    """
    hsp_trajectories_paths = []
    client_trajectories_pahts = []
    #for tp_folder in timepoint_folders:
    #new_input = f'{input_folder}{tp_folder}/'
    
    #filters for JUST colocalised molecules
    yikes = [[f'{root}/{name}' for name in files if '_colocal_traj.csv' in name]for root, dirs, files in os.walk(f'{inputt}/')]
    #flatten list
    yikes = [item for sublist in yikes for item in sublist if tp_folder in item]
    clients = [path for path in yikes if client in path]
    hsp_trajs = [path for path in yikes if hsp in path]
    hsp_trajectories_paths.append(hsp_trajs)
    client_trajectories_pahts.append(clients)

    return hsp_trajectories_paths, client_trajectories_pahts

def fix_backslash(listo):
    l=[]
    for item in listo:
        item=item.replace('\\', '/')
        item=item.replace('//','/')
        l.append(item)
    return l

def read_trajectories(hsp_traj, client_traj):
    """reads in trajectories for clients and hsps, and matches them (i.e. tells you they are in complexes)

    Args:
        hsp_traj (list): this is a list of paths, to files that contain trajectories
        client_traj (list  ): same as hsp, but for client folder

    Returns:
        df: new dataframe with all trajectories, with names indicating whether they are matched or not
    """
    new = []
    for trajectory in hsp_traj:
        #read in trajectory path
        hsp_trajectories = pd.read_csv(trajectory)
        hsp_trajectories.columns = [
            str(col) + '_hsp' for col in hsp_trajectories.columns]
        for traj in client_traj:
            #this line figures out if the two trajectories are from the same TIMEPOINT and if they are, append hsp or client to the end of every molecule name (trajectory name i.e. column name)
            if traj.split('/')[-4] == trajectory.split('/')[-4]:
                client_trajectories = pd.read_csv(traj)
                client_trajectories.columns = [
                    str(col) + '_client' for col in client_trajectories.columns]
                

        for item in hsp_trajectories:
            trajectory_hsp = hsp_trajectories[[item]]
            hspnumber = item.split('_')[1]
            for item1 in client_trajectories:
                client_trajectory = client_trajectories[[item1]]
                clientnumber = item1.split('_')[1]
                #do they match?
                if hspnumber == clientnumber:
                    #if they do match, join them together
                    newtraj = trajectory_hsp.join(client_trajectory)
                    newtraj = newtraj.T
                    newtraj['path_info'] = trajectory
                    #assign a random number so they don't get duplicated
                    newtraj['identifier'] = random.randint(1, 10000000)
                    newtraj['timepoint (min)'] = trajectory.split('/')[-4].split('_')[-1]

                    newtraj = newtraj.reset_index().rename(
                        columns={'index': 'molecule_number'})
                    
                    new.append(newtraj)
    new = pd.concat(new)
    new['identifier']=new['identifier'].astype(str)
    new['new_name'] = new['molecule_number']+'_' + new['identifier']+ '_' +new['timepoint (min)']
    #assign new molecule name
    new['molecule_number'] = [f'{molecule_number}_{x}' for x, molecule_number in enumerate(new['new_name'])]

    return new

def make_directory(o):
    if not os.path.exists(o):
        os.makedirs(o)

def reshape_data(df):

    df = df.T
    df.columns = df.iloc[0]
    df = df.tail(df.shape[0] - 1)
    return df


if __name__ == "__main__":    

    inputt= 'data/example_raw_trajectory_data/'
    output_folder = 'python_results/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    client = 'Client'
    hsp = 'Hsp'

    timepoint_folders = [folder for folder in os.listdir(f'{inputt}')if 'beads' not in folder]

    for tp_folder in timepoint_folders:
        input_folder=f'{inputt}{tp_folder}/'
        #collect all trajectories files from hsp and client
        hsp_traj, client_traj=grab_trajectories_paths_only(inputt, client, hsp, tp_folder)
        hsp_trajs = [item for sublist in hsp_traj for item in sublist]
        client_traj = [item for sublist in client_traj for item in sublist]
        #fix slashes so they match and can be split upon
        hsp_traj=fix_backslash(listo=hsp_trajs)
        client_traj=fix_backslash(listo=client_traj)
        #read, clean, match trajectories, give new names
        named_trajectories=read_trajectories(hsp_traj, client_traj)
        #save this version as a checkpoint  
        named_trajectories.to_csv(f'{output_folder}renamed_trajectories_all.csv')
        #make a list of the columns that ARENT molecule_name and drop all of them!!! because this is the ONLY col you can have in py4bleaching
        test=['identifier','path_info','timepoint (min)','new_name']
        testo=named_trajectories.drop(columns=[col for col in named_trajectories if col in test])

        #now we want to FILTER and resave the trajectories as separate dataframes (remembering they're different fluorophores and as such neeeded to be analysed for photobleaching spearately.)
        #define outputs that nests them correctly for feeding into py4bleaching
        hsp_output = f'{output_folder}/shsp/{tp_folder}/coloc/'
        client_output = f'{output_folder}/client/{tp_folder}/coloc/'
        #create these directories
        make_directory(o=hsp_output)
        make_directory(o=client_output)

        #filter for hsp or client only based on the molecule name, then reshape the data so it fits with the way we read it in in py4bleaching
        hsp_only = testo[~testo['molecule_number'].str.contains('client')]
        hsp_only=reshape_data(df=hsp_only)
        client_only = testo[~testo['molecule_number'].str.contains('hsp')]
        client_only=reshape_data(df=client_only)
        #save them to the outputs you just made
        hsp_only.to_csv(f'{hsp_output}hsp_coloc_traj.csv')
        client_only.to_csv(f'{client_output}client_coloc_traj.csv')

