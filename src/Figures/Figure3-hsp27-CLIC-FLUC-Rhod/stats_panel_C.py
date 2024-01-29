"""script to test how the distribution of molar ratios WITHIN complexes changes through time (this is the data shown in the hexbin plots, specifically Figure 3, panel c i-iii). This uses the Kruskal wallis test to test the distributions, then Dunn's post-hoc test.
"""
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
from scipy import stats
import scikit_posthocs as sp

def kruskal_wallis(all_ratios_df, pair):
    """performs test & post hoc (dunns) on the distribution of RATIOS between sHsps and clients at EACH timepoint in the incubation

    Args:
        all_ratios_df (df): dataframe with all ratios for every complex
        pair (string): the pair of proteins in the complex

    Returns:
        df: a dataframe with the p values for each timepoint (each row is a timepoint)
    """
    #split this up and turn each timepont into an array, so that we can compare between them easily
    zero_tp_group1 = all_ratios_df[all_ratios_df['timepoint (min)']==0.]
    zero_tp_group1 = np.array(zero_tp_group1['ratio'])

    twenty_tp_group1 = all_ratios_df[all_ratios_df['timepoint (min)']==20.]
    twenty_tp_group1 = np.array(twenty_tp_group1['ratio'])

    forty_tp_group1 = all_ratios_df[all_ratios_df['timepoint (min)']==40.]
    forty_tp_group1 = np.array(forty_tp_group1['ratio'])

    sixty_tp_group1 = all_ratios_df[all_ratios_df['timepoint (min)']==60.]
    sixty_tp_group1 = np.array(sixty_tp_group1['ratio'])

    threeh_tp_group1 = all_ratios_df[all_ratios_df['timepoint (min)']==240.]
    threeh_tp_group1 = np.array(threeh_tp_group1['ratio'])

    sevenh_tp_group1 = all_ratios_df[all_ratios_df['timepoint (min)']==420.]
    sevenh_tp_group1 = np.array(sevenh_tp_group1['ratio'])


    comparisons_list = [zero_tp_group1,twenty_tp_group1, forty_tp_group1,sixty_tp_group1, threeh_tp_group1, sevenh_tp_group1]
    #this part tells us whether there is a difference in some of the means of these data sets.
    result = stats.kruskal(zero_tp_group1, twenty_tp_group1, forty_tp_group1, sixty_tp_group1, threeh_tp_group1)
    print(result)

    #if p<0.05 for these comparisons, we need to perform Dunn's test for multiple comparisons on the same datasets
    data = [all_ratios_df[all_ratios_df['timepoint (min)']==0.0]['ratio'], 
            all_ratios_df[all_ratios_df['timepoint (min)']==20.0]['ratio'],
            all_ratios_df[all_ratios_df['timepoint (min)']==40.0]['ratio'],
            all_ratios_df[all_ratios_df['timepoint (min)']==60.0]['ratio'],
            all_ratios_df[all_ratios_df['timepoint (min)']==240.0]['ratio'],
            all_ratios_df[all_ratios_df['timepoint (min)']==420.0]['ratio']]
    p_values = sp.posthoc_dunn(data, p_adjust='holm')

    #now we have a dataframe where every row is a timepoint, and compares every other timepoint to itself. each column the one being compared moves across. To get a simpler output, we can just say p_values < 0.05, and any entry that is statistically significant will be 'True'.
    print(pair)


    for timepoint, df in all_ratios_df.groupby('timepoint (min)'):
        print(timepoint, df['ratio'].median())

    p_values < 0.05
    p_values['pair'] = pair
    return p_values

if __name__ == "__main__":
    #location of the previously saved dataframe which contains your actual 'ratio' columns (stoichiometry at all the timepoints)
    input_folder = ('data/Figures/Figure_3/C-matched_hexbins/hsp27/')
    #make this the same as what it was in the previous script for plotting as this is what it's saved as and we will look for that
    pairs = ['CLIC and hsp27', 'FLUC and hsp27', 'Rhodanese and hsp27']

    st_summary = []
    for pair in pairs:
        all_ratios_df = pd.read_csv(f'{input_folder}{pair}_ratios_added.csv')
        st = kruskal_wallis(all_ratios_df, pair)
        st_summary.append(st)
    st_summary = pd.concat(st_summary)
