"""finds from colocalisation data, the % of client proteins which are colocalised with a sHsp. Plots over time. This script needs to have both NON-COLOC and COLOC molecules run through py4bleaching to be accurate. If not, it will give out 100% colocalisation as the non-colocalised molecules will be excluded.


"""

import pandas as pd
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

def collate(data_files):
    """this function acts to gather the cleaned trajectories that came from py4bleaching, and to split up the molecule name into metadata that is easier to filter for later.

    Args:
        data_files (list): list of files called 'cleaned_data', from two colour experiments, one for each protein (client or hsp)
    """
    molecules = []
    for filepath in data_files:
        filepath = filepath.replace('\\', '/')
        data = pd.read_csv(filepath)
        data.drop([col for col in data.columns.tolist()
                   if ' ' in col], axis=1, inplace=True)

        data[['unique_num', 'timepoint', 'colocalisation', 'protein', 'num']
             ] = data['molecule_number'].str.split('_', expand=True)

        molecules.append(data)

    molecules = pd.concat(molecules)
    return molecules

def calculate_sHsp_coloc(timepoints, shsp, zero, fifteen, thirty, sixty, three_hour, four_hour, seven_hour):
    """for each timepoint find the total number of HSPS (this script is not client molecules) that are colocalised and report as a %. append this to a list for each timepoint

    Args:
        timepoints (list): timepoints
        shsp (str): shsp of interes
        zero (list): these are all lists that can be appended to
        fifteen (_type_): _description_
        thirty (_type_): _description_
        sixty (_type_): _description_
        three_hour (_type_): _description_
        four_hour (_type_): _description_
        seven_hour (_type_): _description_

    Returns:
        list: a bunch of lists with coloc data in them for each tp
    """
    #
    #find % colocalisation for shsps over time
    for timepoint in timepoints:
        timepoint
        protein = 'shsp'
        timepoint_shsp = shsp[shsp["timepoint"] == timepoint]
        total_shsp_timepoint = len(timepoint_shsp)
        coloc_timepoint_shsp = len(
            timepoint_shsp[timepoint_shsp["colocalisation"] == "Coloc"])
        if total_shsp_timepoint > 0:
            percent_colocal_shsp_timepoint = coloc_timepoint_shsp/total_shsp_timepoint*100
        else:
            percent_colocal_shsp_timepoint = 0
        listo = [timepoint, protein, total_shsp_timepoint,
                 coloc_timepoint_shsp, percent_colocal_shsp_timepoint]
        if 'zero' in timepoint:
            zero.append(listo)
        elif '15min' in timepoint:
            fifteen.append(listo)
        elif '30min' in timepoint:
            thirty.append(listo)
        elif '60min' in timepoint:
            sixty.append(listo)
        elif '180min' in timepoint:
            three_hour.append(listo)
        elif '240min' in timepoint:
            four_hour.append(listo)
        elif '420min' in timepoint:
            seven_hour.append(listo)
    return zero, fifteen, thirty, sixty, three_hour, four_hour, seven_hour

def calculate_client_coloc(timepoints, client, zero, fifteen, thirty, sixty, three_hour, four_hour, seven_hour):
    """for each timepoint find the total number of CLIENTS  that are colocalised and report as a %. append this to a list for each timepoint

    Args:
        timepoints (list): timepoints
        shsp (str): shsp of interes
        zero (list): these are all lists that can be appended to
        fifteen (_type_): _description_
        thirty (_type_): _description_
        sixty (_type_): _description_
        three_hour (_type_): _description_
        four_hour (_type_): _description_
        seven_hour (_type_): _description_

    Returns:
        list: a bunch of lists with coloc data in them for each tp
    """
    #find % colocalisation for client over time
    for timepoint in timepoints:
        timepoint
        protein = 'client'
        timepoint_client = client[client["timepoint"] == timepoint]
        total_client_timepoint = len(timepoint_client)
        coloc_timepoint_client = len(
            timepoint_client[timepoint_client["colocalisation"] == "Coloc"])
        percent_colocal_client_timepoint = coloc_timepoint_client/total_client_timepoint*100
        listo = [timepoint, protein, total_client_timepoint,
                 coloc_timepoint_client, percent_colocal_client_timepoint]
        if 'zero' in timepoint:
            zero.append(listo)
        elif '15min' in timepoint:
            fifteen.append(listo)
        elif '30min' in timepoint:
            thirty.append(listo)
        elif '60min' in timepoint:
            sixty.append(listo)
        elif '240min' in timepoint:
            four_hour.append(listo)
        elif '420min' in timepoint:
            seven_hour.append(listo)
    return zero, fifteen, thirty, sixty, three_hour, four_hour, seven_hour

def scatter_plot_client(for_plotting, output_folder):
    """plot CLIENT colocalisation at each timepoint as a lineplot

    Args:
        for_plotting (df): df with the colocalisation
        output_folder (str): where to save
    """
    #filter for only the client colocalisation
    client_only = for_plotting[for_plotting['protein'] == 'client']
    #plot as a scatterplot
    ax = plt.scatter(x=client_only['timepoint'],
                     y=client_only['value'],
                     c='darkviolet')
    ax = plt.plot(client_only['timepoint'], client_only['value'], c='plum')
    ax = plt.gca()
    ax.set_ylim([0, 100])
    plt.ylabel('Percent colocalisation (%)')
    plt.title('Percentage of client colocalised with sHsp')
    plt.xlabel('Time (h)')
    ax = ax.get_figure()
    ax.savefig(f'{output_folder}scatter_percent_colocal')


if __name__ == "__main__":    
    input_folder = 'data/'
    output_folder = 'python_results/colocalisation/'

    if not os.path.exists(output_folder):
         os.makedirs(output_folder)

    #this dataframe contains ALL trajectories from both proteins, and the metadata so I can grab specific ones later if required
    #find the df containing cleaned colocalisation data
    data_files = [[f'{root}/{filename}' for filename in files if 'cleaned_data.csv' in filename]for root, dirs, files in os.walk(f'{input_folder}')]

    data_files = [item for sublist in data_files for item in sublist if 'two-colour' in item]

    molecules = collate(data_files)
    timepoints = molecules['timepoint'].unique().tolist()


    #split big df into hsp and client again
    shsp = molecules[molecules["protein"]=="hsp"]
    client = molecules[molecules['protein']=='client']
    zero = []
    fifteen = []
    thirty = []
    sixty = []
    three_hour = []
    four_hour = []
    seven_hour = []


    zero, fifteen, thirty, sixty, three_hour, four_hour, seven_hour = calculate_sHsp_coloc(
        timepoints, shsp, zero, fifteen, thirty, sixty, three_hour, four_hour, seven_hour)


    zero, fifteen, thirty, sixty, three_hour, four_hour, seven_hour = calculate_client_coloc(
        timepoints, client, zero, fifteen, thirty, sixty, three_hour, four_hour, seven_hour)

    percent_colocalisation_all = [zero, fifteen,
                                thirty, sixty, four_hour, seven_hour]
    #list names of columns
    column_names = ['timepoint', 'protein', 'total number of proteins',
                    'colocalised proteins', 'percent colocalisation']
    percent_colocalisation_all = pd.DataFrame([item for sublist in percent_colocalisation_all for item in sublist])
    #assign column names based on list defined
    percent_colocalisation_all.columns = column_names
    #save it
    percent_colocalisation_all.to_csv(f'{output_folder}percent_colocalisation_all.csv')

    #changes the timepoint column from minutes into integers that represent hours, so that the scatter plot is scaled for hours on the x axis correctly (rather than equal distance between timepoints)
    timepoint_dict = {'zero': 0, '15min': 0.25,
                    '30min': 0.5, '240min': 4, '60min': 1, '420min': 7}
    percent_colocalisation_all['timepoint'] = percent_colocalisation_all['timepoint'].map(timepoint_dict)
    for_plotting = pd.melt(percent_colocalisation_all, id_vars=['timepoint', 'protein'], value_vars='percent_colocalisation')

    #this plots the colocalisation for this experiment
    scatter_plot_client(for_plotting, output_folder)
