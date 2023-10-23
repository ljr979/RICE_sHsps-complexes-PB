from src.utilities.pb.py4bleaching import analysis

#input is the path to the imagejresults folder containing trajectories files. output_folder is python_results and whatever treatment you want it to be nested within
input_folder = 'data/example_python_output/two-colour/0_matching_trajectories/'
output_folder = 'python_results/py4bleaching/'

model_name = 'Model_2'


analysis.pipeline(input_folder, output_folder, probability_threshold=0.7, model_name=model_name, x_norm=False, matched=True)

