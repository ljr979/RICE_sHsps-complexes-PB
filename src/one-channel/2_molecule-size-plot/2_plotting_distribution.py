
#below script is to PLOT the distributions of the molecule sizes in a violin plot
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
#input folder is the folder that the py4bleaching just saved
input_folder = 'python_results/'

output_folder = 'python_results/'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

#gather the molecule count files for all the time points
stoich_files = [[f'{root}/{filename}' for filename in files if 'molecule_counts.csv' in filename]for root, dirs, files in os.walk(f'{input_folder}')]
#flatten list
stoich_files = [item for sublist in stoich_files for item in sublist if 'one-colour' in item]

#concatinate any of the same name (molecule_counts) into one big file. This would be e.g. from different time points that are in separate files etc.
def concatinate_data(stoich_files):

    molecule_counts = []

    for filepath in stoich_files:
        data = pd.read_csv(filepath)
        tp = data.treatment
        timepoints = []
        for item in tp:
            timepoint = item.split('-')[0]
            timepoints.append(timepoint)
        data['timepoint'] = timepoints

        molecule_counts.append(data)
    molecule_counts = pd.concat(molecule_counts)

    drop_cols = ['all_small_mol_count','single_step_mol_count', 'max_fluorescence']
    #we also drop columns here that are not needed for this analysis (we use last step fluorescence and as such drop the others)
    molecule_counts.drop([col for col in molecule_counts.columns.tolist(
    ) if ' ' in col], axis=1, inplace=True)
    molecule_counts.drop([col for col in molecule_counts.columns.tolist(
    ) if col in drop_cols], axis=1, inplace=True)
    return molecule_counts

molecule_counts=concatinate_data(stoich_files)
df = pd.melt(molecule_counts, id_vars=['protein', 'colocalisation', 'molecule_number',
                                       'treatment', 'timepoint'], value_vars=['last_step_mol_count'])

ax = sns.violinplot(x="timepoint", y="value",  data=df, order=[
                    'zero', '20min', '40min', '60min'], scale='width', palette='viridis')

ax.set_ylabel('# of subunits')
ax.set_ylim(0, max(df['last_step_mol_count']))
plt.title(f'# of subunits/fluorescent spot')
plt.savefig(f'{output_folder}_stoichiometries.png')
plt.show()

