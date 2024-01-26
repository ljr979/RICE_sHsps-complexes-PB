
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
from scipy import stats
import scikit_posthocs as sp
#now do same stats on the start and finish data for each client, but just kruskal wallis
#ran in the data from aBc_Fig3_plotting up to the point where I save them ,so could alternatively just read in here the dfs 'Clic.csv' and 'fluc.csv' and aB-c.csv
def kruskal_wallis_start_finish(treatment, df):
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

source_data = 'python_results/aBc_Fig3_plotting/'
clients=['CLIC', 'FLUC', 'RHOD']

proteins=['CLIC', 'FLUC', 'aBc']

for x in proteins:
    protein=x
    #for CLIC, read in df
    df = pd.read_csv(f'{source_data}{x}.csv')
    #lets first look at the data where CLIC was incubated with aB-c
    if x in clients:
        treatment = '+sHsp'
        #now we are going to store the variables that are the arrays of the start and end for this treatment, so that we can compare further later on
        print(x)
        start_hsp, end_hsp =kruskal_wallis_start_finish(treatment, df=df)

        treatment = '-sHsp'

        #now we are going to store the variables that are the arrays of the start and end for this treatment, so that we can compare further later on
        start_nohsp, end_nohsp =kruskal_wallis_start_finish(treatment, df=df)

        #now we compare between the end of the incubation for the two treatments (did the size of the client vary after 7h when we had sHsp present or not?)
        result_end_each = stats.kruskal(end_hsp, end_nohsp)
        print('comparing END timepoint mol (client)size distribution between nosHsp present, and a shsp present:')
        print(result_end_each)
    if x not in clients:
        treatment = '-client'
        #now we are going to store the variables that are the arrays of the start and end for this treatment, so that we can compare further later on
        print(x)
        start_noclient, end_noclient = kruskal_wallis_start_finish(treatment, df=df)

        treatment = '+client'

        pairs=df[df['incubated']==treatment]['Pair'].unique().tolist()

        for client in pairs:
            client
            print(client)
            df1=df[df['Pair']==client] 
            #now we are going to store the variables that are the arrays of the start and end for this treatment, so that we can compare further later on
            start_client, end_client = kruskal_wallis_start_finish(treatment, df=df1)


            #now we compare between the end of the incubation for the two treatments (did the size of the client vary after 7h when we had sHsp present or not?)
            result_end_each = stats.kruskal(end_noclient, end_client)
            print(f'comparing END timepoint mol (sHsp) size distribution between no {client} and a {client} present:')
            print(result_end_each)






