"""this script performs two-way anova on the molecule size data in Figure 4 (panels c-d), and Supp Figure 4 (panels a & c). This compares between colocalised and non-colocalised molecules, and between all timepoints. These were then annotated in the figures.

Returns:
    df: dataframe with all the p values for each protein, timepoint, and colocalisation state
"""

import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
from scipy import stats
import scikit_posthocs as sp
import statsmodels.stats.multicomp as mc
import statsmodels.api as sm
from statsmodels.formula.api import ols
from bioinfokit.analys import stat

def two_way_ANOVA(d, pval, variable='last_step_mol_count', xvar1='Timepoint', xvar2='Colocalisation'):
    """performs two-way ANOVA on data, with Tukey's post-hoc

    Args:
        d (df): dataframe with all molecule size data at all timepoints
        pval (float): p value deemed significant
        variable (str, optional): the variable you want to compare. Defaults to 'last_step_mol_count'.
        xvar1 (str, optional): important variable 1 . Defaults to 'Timepoint'.
        xvar2 (str, optional): important variable 2. Defaults to 'Colocalisation'.

    Returns:
        df: dataframe with all of the p values for each tp and coloc state.
    """
    two_sigs = []
    for protein, df in d.groupby('Protein'):
        protein 
        df
        print(protein)
        model = ols(f'{variable} ~ C({xvar1}) + C({xvar2}) + C({xvar1}):C({xvar2})', data=df).fit()
        anova_table = sm.stats.anova_lm(model, typ=3)
        anova_table
        # perform multiple pairwise comparison (Tukey HSD)
        # unequal sample size data, tukey_hsd uses Tukey-Kramer test
        res = stat()
        # for main effect 
        res.tukey_hsd(df=df, res_var=variable, xfac_var=xvar1, anova_model=f'{variable}~C({xvar1})+C({xvar2})+C({xvar1}):C({xvar2})')
        x = res.tukey_summary
        sigx = x[x['p-value']<=pval]
        res.tukey_hsd(df=df, res_var=variable, xfac_var=[f'{xvar1}',f'{xvar2}'], anova_model=f'{variable}~C({xvar1})+C({xvar2})+C({xvar1}):C({xvar2})')
        y = res.tukey_summary
        ysig = y[y['p-value']<=pval]
        ysig['Protein'] = protein
        two_sigs.append(ysig)

    two_sigs = pd.concat(two_sigs)
    return two_sigs

if __name__ == "__main__":

    output_folder = 'data/Figures/Figure_3/aBc/'
    data = pd.read_csv('data/Figures/Figure_3/aBc/aBc_coloc_noncoloc_alltps.csv')
    proteins = ['FLUC', 'CLIC']

    for b in proteins:     
        d = data[data['Pair']==f'{b}_aB-c']
        two_sigs_b = two_way_ANOVA(d=d, pval=0.05, variable='last_step_mol_count', xvar1='Timepoint', xvar2='Colocalisation')
        two_sigs_b.to_csv(f'{output_folder}twoway_ANOVA_{b}.csv')