# Data directory
```Images/```

**complexes**

*additional-for-imagej*
- This folder contains the files that need to exist for the imageJ pre-processing to work. This is a Gaussian blur of each channel to subtract during beam profile correction, an ROI file to define the matched channels, and the transformation matrices to align the channels for colocalisation
  
*Client-AF647*
- Example image of the 647nm emission obtained from FLUC-AF647 incubated with Hsp27


*Hsp-AF488*
- Matching 488nm emission from Hsp27 AF488. Used to find colocalisation (complexes) between the two proteins

-----------------
**controls**

*Client-AF647/*

*Hsp-AF488/*

- Example image within each of these folders of FLUC (647 nm) or Hsp27 (488 nm) incubated alone and photobleached
------------------------------------------------
```Imagej_data/```

**complexes**
*0min*, *20min*

- This section pertains to **two-colour** experiments only (where both colours have been imaged in the same experiment, for complex formation)
- Within each timepoint folder is ```non-coloc``` and ```coloc``` folders, containing subfolders ```Client``` and ```Hsp```. Each of these folders has an example ImageJ trajectories output for either colocalised molecules or non-colocalsied molecules. Each 'traj' file contains the fluorescence trajectory time traces for every molecule in an image (i.e. one file = one image)
- These folders should act as an example of the appropriate folder nesting to feed trajectories data to py4bleaching. 

**controls**
*Client-AF647*, *Hsp-AF488*

- Each of these folders contains the output from the ImageJ scripts pre-processing the controls images (single colour, **either** a client or a chaperone imaged, not together).
- This includes a file containing the x and y position of each molecule detected, and a file with the fluorescence trajectory of each molecule, to be fed to py4bleaching.

---

```example_python_output/```

**one-colour**
- This section is an example of the output we would obtain after running the 'py4bleaching' utilities package on an experiment with one fluorescently labelled protein. 
  
- Specifically, in the 'one-colour'section the file (FLUC_molecule_counts.csv) contains the number of the AF647-labelled FLUC molecules per fluorescent spot, at time points ('treatment') over the course of an hour. This is gained from the type of example data in the 'imagej_data/controls/Client-AF647/' folder, in the trajectory files, after they've been processed by py4bleaching.

- This type of data can be plotted as a violinplot to check for distributions of molecule sizes using the script src/one-channel/molecule-size-plot/2_plotting_distribution.py. 


**two-colour**

*0_matching_trajectories*
- This section contains example output from the script ```0_match_trajectories.py``` (in ```src/two-channels/1_trajectory-analysis/```)
- This comes from the ```imagej_data/complexes/``` data section, which has two timepoints included. As such, the example python output is shown for both timepoints and proteins.
- These files (i.e. ```client_coloc_traj.csv``` and ```hsp_coloc_traj.csv```) contain the original trajectories, but have been renamed such that their 'molecule name' contains a unique number, some metadata about the timepoint, and this information and number matches the corresponding sHsp or client molecules that is within the complex. 
  
- These are nested appropriately within folders for py4bleaching analysis.

*1_trajectory_analysis*
- This section contains trajectories for a sHsp and a client protein which were incubated together to form complexes with their names matched as described in ```0_matching_trajectories``` (these are the same as those in the above folder, but are placed here for ease of use also)
- py4bleaching is run separately on the raw trajectories from different fluorophores (i.e. don't run it on sHsp and client at the same time, rather direct py4bleaching to e.g. ```1_trajectory_analysis/client/```). Importantly, this should ALSO have a non-colocalised molecules folder with the trajectories of the non-colocalised molecules. The data provided is an example only, but in reality, there would be a subfolder under ```client/``` and ```hsp/``` , and each timepoint, which was ```non-coloc``` (as in the folder ```coloc```). This is necessary if you are going to run the colocalisation script on the py4bleaching output. 
- These files are used to do colocalisation analysis on (```src/two-channels/1_colocalisation```).
- Under each protein subheading there is a folder called 'py4bleaching' which has the example output of py4bleaching (only included the full output for client, and 'cleaned_data' for hsp, as this allows colocalisation analysis example to run).
- The ```molecule_counts.csv``` is the number of subunits per molecule (trajectory) throughout this work. 