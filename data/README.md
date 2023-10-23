# Data directory
*Data/*
/example_python_output
*/one-colour*
- this section is an example of the output we would obtain after running the 'py4bleaching' utilities package on an experiment with one fluorescently labelled protein. 
  
- Specifically, in the 'one-colour'section the file (FLUC_molecule_counts.csv) contains the number of the AF647-labelled FLUC molecules per fluorescent spot, at time points ('treatment') over the course of an hour. This is gained from the type of example data in the 'imagej_data/controls/Client-AF647/' folder, in the trajectory files, after they've been processed by py4bleaching.

-This type of data can be plotted as a violinplot to check for distributions of molecule sizes using the script src/one-channel/molecule-size-plot/2_plotting_distribution.py. 


*/two-colour*

0_matching_trajectories/
-this section contains example output from the script 0_match_trajectories.py (in src/two-channels/1_trajectory-analysis/)
-these files (i.e. client_coloc_traj.csv and hsp_coloc_traj.csv) contain the original trajectories, but have been renamed such that their 'molecule name' contains a unique number, some metadata about the timepoint, and this information and number matches the corresponding sHsp or client molecules that is within the complex. 

1_trajectory_analysis/
-This section contains example files 'cleaned_data.csv' for a sHsp and a client protein which were incubated together to form complexes. 
- py4bleaching was run separately on them (as it currently should be run separately when the fluorophores are different), thus they are different files from molecules in the same experiment.
- These files are used to do colocalisation analysis on (src/two-channels/1_colocalisation).
  
----------------------------------------------------------------------------------

*Imagej_data/*

complexes/
-------------------
controls/



----------------------------------------------------------------------------------
*Images/*
complexes/
/additional-for-imagej/
Client-AF647/
Hsp-AF488/
-----------------
controls/
Client-AF647/
Hsp-AF488/








----------------------------------------------------------------------------------

CODE BASE
*src/*
0_imagej-processing/
- this folder contains all of the imagej macros to be run prior to any python code (for image processing). These will be different depending on the experimental design and are split up according to our experimental design (i.e. 488 labelled Hsp alone, one colour experiment(/1_Hsp-only/488_trajectories.ijm) and the example outputs are provided in the data/imagej_data/ controls/ section.)
  
- essentially these macros subtract background, identify molecules, and extract their change in fluorescence over time. Also outputs the X & Y positions of each molecule.
  
- AF647 labelled proteins alone (not incubated with a 488-labelled protein) is the same as above, but in 1_Client-only/647_trajectories.ijm. Outputs in imagej_data/controls/Client-AF647/ section. 

2_Complexes/
- scripts same as 1_Controls macros, but this macro aligns the two channels as well, finds colocalised spots, and then extracts their trajectory over time. 
  
macros-required/
- this folder contains macros and plugins required to execute the macros successfully.
- copy and paste into plugins folder in imagej location on desktop.

one-channel/
1_trajectory_analysis.py
- this script is run to feed your x_traj.csv files to py4bleaching. Send it to the location of the single condition e.g. data/imagej_data/controls/Client-AF647/ AND IT WILL FIN