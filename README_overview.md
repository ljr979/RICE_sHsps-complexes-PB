
# RICE_sHsps-complexes-PB

This repository contains the analysis code associated with the sHsp complexes project, led by Professor Heath Ecroyd. This manuscript has been submitted for publication under the title **"#####"**.

This manuscript has been submitted as a preprint via BioRxiv [here](biorxiv/link). A link to the final version will be provided upon publication.

## Prerequisites

This analysis assumes a standard installation of Python 3 (=> **3.7.10**), and ImageJ =>**v1.53c**. For specific package requirements, see the environment.yml file, or  create a new conda environment containing all packages by running ```conda create -f environment.yml```. In addition to the analysis contained here, the simple statistical tests & data fitting in Figure 1 were performed using [GraphPad Prism v **9.0**](https://www.graphpad.com/scientific-software/prism/).

For ImageJ processing, the required plugins and macros can be found in 0_imagej-processing/macros-required/ and should be transferred into the ImageJ plugins folder.

## Raw data

For convenience, example datasets are provided here under the ```data``` folder. These data may be used to explore the workflows presented here as described below.

In addition, data used to generate each of the figures in this work, have been uploaded as an open-access [Zenodo dataset](https://doi.org/###/zenodo.###). These datasets can be collected automatically using the ```raw_data.py``` script in each of the respective analysis folders.

## Workflow

Initial preprocessing of the raw TIFF images was performed in ImageJ, using the appropriate macros, which can be found in ```src/0_imagej-procesing```. This folder contains all of the imagej macros to be run prior to any python code (for image processing). These will be different depending on the experimental design and are split up according to our experimental design (i.e. 488 labelled Hsp alone, one colour experiment(```/1_Hsp-only/488_trajectories.ijm```) and the example outputs are provided in the ```data/imagej_data/controls/```section). 

Individual analyses are presented within the ```src``` folder, under either ```one-channel``` or ```two-channels``` depending on the nature of the data to be analysed. Where processing order is important for individual analyses, scripts have been numbered and should be run in order before unnumbered counterparts.

In addition, each figure can be generated using the scripts provided under the ```src/Figures/``` folder, and the statistics pertaining to the data in each figure can be reproduced using the scripts within the same folder.

1. [one-channel](src/one-channel)

| Script      | Language/Interpreter | Description   |
|-------------|----------------------|---------------|
| ```647_trajectories.ijm``` or ```488_trajectories.ijm``` | ImageJ        | Subtract background, identify molecules, and extract their fluorescence over time. Also outputs the X & Y positions of each molecule. |
| ```1_trajectory_analysis.py``` | Python               | Gather all trajectories from imageJ output, and find the number of subunits in each molecule from photobleaching step size (see ```py4bleaching``` README for more details) |
| ```2_plotting_distribution.py``` |Python |  Visualise calculated molecule sizes in violinplots|



2. [two-channels](src/two-channels)

| Script      | Language/Interpreter | Description   |
|-------------|----------------------|---------------|
| ```Colocalised_client_hsp-alignment.ijm ```| ImageJ        | Same as ```1_Controls``` macros, but this macro aligns the two channels as well, finds colocalised spots, and then extracts their trajectory over time.  |
| ```0_match_trajectories.py``` | Python               | Finds trajectories from ImageJ and matches those that are colocalised (within complexes) to one another, assigning them a new name so they can be identified as being specifically in complex with each other prior to molecule size analysis |
|```1_trajectory_analysis.py```|Python|Use ```py4bleaching``` to calculate molecule size for all trajectories. Can either run this on the separate files as output from imageJ, or can run on the 'matched trajectories' file from ```0_match_trajectories.py```, using the 'matched' keyword|
|```2_colocalisation_analysis.py```|Python|Find % of client molecules in complexes|
|```3_matched_counts_for_hexbins```|Python|If ```0_match_trajectories.py``` was not run prior to py4bleaching, this script will match trajectories and find their molecule size from the step size analysis, then save them for plotting|
|```3_stoichiometries_kdeplots.py```|Python|Plot complexes as kdeplots to visualise complex size distribution over time|


## References

[1]: my/really/cool/link

1. Reference Paper (2020) Really Cool Paper. **High Impact Journal**. https://doi.org/####/####
