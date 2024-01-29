"""This script statistically tests the difference in the distribution between the molecule size of the clients in the presence and absence of the sHsp (aB-c). 

It tests and then prints the results of this test, which were then annotated on Figure 2 panel C, and Supplementary figure 2. 
"""
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
from scipy import stats
import scikit_posthocs as sp

def kruskal_wallis_start_finish(treatment, df):
    """performs kruskal wallis statistical test on the molecules size distribution, between the start and end of the incubation.

    Args:
        treatment (str): string which says which treatment you want to compare the start and end between.
        df (df): dataframe to subset on the treatment

    Returns:
        array: two arrays, which have the start and end data. also PRINTS the results of the test. 
    """
    #this function performs the kruskal wallis test on the start and finish of incubation of a particular protein (could be client or sHsp). it filters the original df for the treatment you want (+ or - client or sHsp), then it grabs the molecule counts of the start and the end for that treatment. 
    #then, it does the test and returns the result as well as the start and end arrays.
    #so df should be the protein.csv file saved in the repo
    #treatment should be a string, plus or minus client or shsp
    print(f'{treatment}:')
    treatment=df[df['incubated'] == treatment]
    start = np.array(treatment[treatment['start_end'] == 'start']['last_step_mol_count'])
    end = np.array(treatment[treatment['start_end'] == 'end']['last_step_mol_count'])
    results = stats.kruskal(start, end)
    
    print(results)

    return start, end
if __name__ == "__main__":
    # read in here the dfs 'Clic.csv' and 'fluc.csv' and aB-c.csv (these have the number of molecules for each of these proteins)
    input_folder='data/Figures/Figure_2/C-violinplots/'
    clients=['CLIC', 'FLUC', 'RHOD']
    proteins=['CLIC', 'FLUC', 'aBc']

    for x in proteins:
        protein=x
        # read in df for client
        df = pd.read_csv(f'{input_folder}{x}.csv')
        #lets first look at the data where CLIC was incubated with aB-c
        if x in clients:
            #now we are going to store the variables that are the arrays of the start and end for this treatment, so that we can compare further later on
            print(x)
            start_hsp, end_hsp=kruskal_wallis_start_finish(treatment='+sHsp', df=df)

            #now we are going to store the variables that are the arrays of the start and end for this treatment, so that we can compare further later on
            start_nohsp, end_nohsp =kruskal_wallis_start_finish(treatment='-sHsp', df=df)
            #now we compare between the end of the incubation for the two treatments (did the size of the client vary after 7h when we had sHsp present or not?)
            result_end_each = stats.kruskal(end_hsp, end_nohsp)
            print('comparing END timepoint mol (client)size distribution between nosHsp present, and a shsp present:')
            print(result_end_each)
        if x not in clients:
            treatment = '-client'
            #now we are going to store the variables that are the arrays of the start and end for this treatment, so that we can compare further later on
            print(x)
            start_noclient, end_noclient = kruskal_wallis_start_finish(treatment='-client', df=df)
            pairs=df[df['incubated']=='+client']['Pair'].unique().tolist()

            for client in pairs:
                client
                print(client)
                df1=df[df['Pair']==client] 
                #now we are going to store the variables that are the arrays of the start and end for this treatment, so that we can compare further later on
                start_client, end_client = kruskal_wallis_start_finish(treatment='+client', df=df1)
                #now we compare between the end of the incubation for the two treatments (did the size of the client vary after 7h when we had sHsp present or not?)
                result_end_each = stats.kruskal(end_noclient, end_client)
                print(f'comparing END timepoint mol (sHsp) size distribution between no {client} and a {client} present:')
                print(result_end_each)






