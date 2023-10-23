macro "traj macro [h]" {
	
//this macro is suited to complexes data gathered with separate imagining conditions. 
//that is, if molecules are FRET ing, and you have to bleach them separately, but then look for colocalisation
//this is opposed to dual imaging 
//hsp setup- define the image to get data from
//select hsp first, then client next.
var abc_path = File.openDialog("RAW image hsp");
var dir = File.getParent(abc_path);
var name = File.getName(abc_path);

splitter = split(dir, File.separator)
Array.print(splitter)
len = splitter.length
trim = len -1
ROIpatharr = Array.slice(splitter,0,trim)
Array.print(ROIpatharr)
str = "";
filesep = File.separator;
     for (i=0; i<ROIpatharr.length-1; i++) 
         str = str + ROIpatharr[i] + filesep; 
     str = str + ROIpatharr[ROIpatharr.length-1] + filesep; 
ROIpath = filesep + str + "ROI.zip";
print("abc Raw Image path= " + abc_path)
print("ROI path= " + ROIpath);
output = dir + filesep
print("Output path= " + output)
transformationfile = filesep + str + "TransformationMatrices.txt";

run("Set Measurements...", "centroid stack redirect=None decimal=3");


//client setup

var client_path = File.openDialog("RAW image client");
var dir = File.getParent(client_path);
var name = File.getName(client_path);

splitter = split(dir, File.separator)
Array.print(splitter)
len = splitter.length
trim = len -1
ROIpatharr = Array.slice(splitter,0,trim)
Array.print(ROIpatharr)
str = "";
filesep = File.separator;
     for (i=0; i<ROIpatharr.length-1; i++) 
         str = str + ROIpatharr[i] + filesep; 
     str = str + ROIpatharr[ROIpatharr.length-1] + filesep; 
ROIpath = filesep + str + "ROI.zip";
print("Raw Image path= " + client_path)
print("ROI path= " + ROIpath);
output = dir + filesep
print("Output path= " + output)
//this is the file we generate to adjust for any misalignment between the emission channels.generated using optical beads image
transformationfile = filesep + str + "TransformationMatrices.txt";

//create output folder for all Trajectories (top folder, then beneath this and at each section below, will add in Hsp vs Client, 
//and colocalisation file
trajectories_output = filesep + str + "/Trajectories/"
print(trajectories_output);
File.makeDirectory(trajectories_output);

run("Set Measurements...", "centroid stack redirect=None decimal=3");


client_results = output + "client_results.csv" // this is the path for where the table results are
Hsp_results = output + "Hsp_results.csv"
client_colocalisation_results = output + "client_colocalisation.csv"

Hsp_MAX = output + "MAX_Hsp peaks.tif"
Hsp_image = output + "Hsp.tif"
client_image = output + "client background corrected.tif"

Hsp_background_image = filesep + str + "488_background.tif"
client_background_image = filesep + str + "647_background.tif"

//now that everything is defined we can start doing the processing.
open(client_path);	//open original raw image
Imagetitle =getTitle()
rename("Raw");

roiManager("Open", ROIpath);
roiManager("Select", 0);	//this will be the client channel
run("Duplicate...", "title=client duplicate");
	

//Background correction 
open(client_background_image);
rename("AVG_client");
selectWindow("client");
run("32-bit");
run("Beam Profile Correction", "electronic_offset=1500 background_image=AVG_client");
setMinAndMax(0, 65536);
run("16-bit");

//crop second time in order to remove black edges in channel
roiManager("Select", 2);
run("Duplicate...", "title=[client background corrected] duplicate");
// this will change the size of the pixels
saveAs("tiff", output + "client background corrected.tif");	//this will save the movie that can be used to get the intensity trajectories
//z stack max projection to get an image where we can get coordinates of all molecules for later extracting their trajectory
run("Z Project...", "stop=20 projection=[Max Intensity]");
saveAs("tiff", output + "MAX_client.tif");
roiManager("reset");
selectWindow("MAX_client.tif");
run("Enhance Contrast", "saturated=0.35");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");

//Find the peaks within the client image
run("Peak Finder");// need to check the threshold here
waitForUser("Have you selected your peaks? \n \nPress OK to continue");
roiManager("Measure");	//this will allow for a table to be produced
saveAs("Results", output + "client_results.csv"); 
close();
selectWindow("AVG_client");
saveAs("jpeg", output + "AVG_client.jpeg");
close();
selectWindow("client background corrected.tif");
close()
selectWindow("client");
close()
selectWindow("Results"); 
     run("Close"); 
roiManager("reset");




//now onto the Hsp channel- this is a little bit tricker as you need to align the channel using bead matrix
open(abc_path);	//open original raw image
rename("Raw hsp");
selectWindow("Raw hsp")

roiManager("Open", ROIpath);
roiManager("Select", 1);	//this will be the hsp (488-labelled) channel
run("Duplicate...", "title=Hsp duplicate");

//Background correction 
open(Hsp_background_image);
rename("AVG_Hsp");
selectWindow("Hsp");
run("32-bit");
run("Beam Profile Correction", "electronic_offset=700 background_image=AVG_Hsp");
setMinAndMax(0, 65536);
run("16-bit");

selectWindow("AVG_Hsp");
saveAs("jpeg", output + "AVG_Hsp.jpeg");
close();
selectWindow("Raw hsp");
close();

//alignment of Hsp channel using bead image
//makes a new directory where the original stack is split into individual tifs ( one for every single frame)
//every image (frame) is aligned to the 647 channel, based on the alignment from the beads
//these aligned images are then saved into the Hsp folder
splitDir= dir + "/Hsp/"
print(splitDir);
File.makeDirectory(splitDir);
list = getFileList(dir); 
selectWindow("Hsp");
run("Image Sequence... ", "format=TIFF use save=splitDir");
selectWindow("Hsp");
close();

outputFolder = splitDir

function action(outputFolder, outputFolder, filename) {
        open(outputFolder + filename);  //+ filesep
        run("MultiStackReg", "stack_1="+filename+" action_1=[Load Transformation File] file_1=["+ transformationfile +"] stack_2=None action_2=Ignore file_2=[] transformation=Affine");
        roiManager("Select", 2);
        run("Duplicate...", "duplicate");
        saveAs("tif", outputFolder + filesep + filename);
        close();
        close();
}


input = "outputFolder";
output = "alignedoutputFolder";

list = getFileList(outputFolder);
for (i = 0; i < list.length; i++)
        action(outputFolder, outputFolder, list[i]);

output = dir + filesep
//now the Hsp folder gets turned into a new, aligned stack, and z projected as we did before.
run("Image Sequence...", "open=[outputFolder] sort");
rename("Hsp");
run("Z Project...", "projection=[Max Intensity]");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");
selectWindow("MAX_Hsp");
saveAs("tiff", output+getTitle()); 
close()
selectWindow("Hsp");
saveAs("tiff", output+getTitle()); 
run("Z Project...", "stop=20 projection=[Max Intensity]");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");
saveAs("tiff", output+getTitle()); 
selectWindow("MAX_Hsp.tif");
run("Enhance Contrast", "saturated=0.35");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");
roiManager("reset");

//Find the peaks within the Hsp image
run("Peak Finder");// need to check the threshold here
waitForUser("Have you selected your peaks? \n \nPress OK to continue");
roiManager("Measure");	//this will allow for a table to be produced
saveAs("Results", output + "Hsp_results.csv"); 
close();
selectWindow("Hsp.tif");
close()
selectWindow("Results"); 
     run("Close"); 
roiManager("reset");

//create colocalisation folders within trajectories folder, for client and hsp

colocalised_output_folder = trajectories_output + "/Coloc/"
print(colocalised_output_folder);
File.makeDirectory(colocalised_output_folder);

//make subfolders for hsp vs client
client_coloc_trajectories = colocalised_output_folder + "/Client/"
print(client_coloc_trajectories);
File.makeDirectory(client_coloc_trajectories);

hsp_coloc_trajectories = colocalised_output_folder + "/Hsp/"
print(hsp_coloc_trajectories);
File.makeDirectory(hsp_coloc_trajectories);

// colocalisation analysis- we look for colocalised spots between the channels (complexes)
//we also save which spots are colocalised and which are not, so we can filter and examine them later
open(client_results);	//this needs to be changed, need to put the path into the begginning of the macro
//selectWindow("client_results.csv");
IJ.renameResults("client");
open(Hsp_results);	//this needs to be changed, need to put the path into the begginning of the macro
//selectWindow("Hsp_results.csv");
IJ.renameResults("Hsp");
run("count colocalized peaks Andrew", "table_1=client table_2=Hsp maximum_distance=3");
selectWindow("client");
IJ.renameResults("Results");
selectWindow("Hsp");
		run("Close")
		// this needs to be the local path of this macro
runMacro(//INSERT PATH HERE: "/fiji-win32/Fiji.app/macros/Filter_for_colocal.ijm")
selectWindow("Colocalization");
saveAs("tiff", output + "Colocalization.tiff");
close()
IJ.renameResults("client_colocalisation")
selectWindow("client_colocalisation")
saveAs("Results", output + "client_colocalisation.csv"); 
selectWindow("client_colocalisation.csv");
     run("Close"); 


open(client_image);
open("client_colocalisation.csv");
IJ.renameResults("Results")
//runMacro("results_to_ROI.txt"); // you need to put the location of your script
runMacro( //INSERT PATH HERE:"/fiji-win32/Fiji.app/macros/results_to_ROI.txt");
selectWindow("Results");
        run("Close");
runMacro(//INSERT PATH HERE:"/fiji-win32/Fiji.app/macros/peak_intensity.txt"); // you need to put the location of your script
selectWindow("Results");
        saveAs("Results", output + "client_colocal_traj.csv");
selectWindow("Results");   
		saveAs("Results", client_coloc_trajectories + Imagetitle +"_Client_colocal_traj.csv");  
        run("Close");
selectWindow("client background corrected.tif");
        run("Close");
roiManager("reset");

open(Hsp_image);	//open original raw image
open("client_colocalisation.csv");
IJ.renameResults("Results")
runMacro(//INSERT PATH HERE:"/fiji-win32/Fiji.app/macros/results_to_ROI_X2.txt"); // you need to put the location of your script
selectWindow("Results");
        run("Close");
runMacro(//INSERT PATH HERE:"/fiji-win32/Fiji.app/macros/peak_intensity.txt"); // you need to put the location of your script
selectWindow("Results");
        saveAs("Results", output + "Hsp_colocal_traj.csv");
selectWindow("Results");
     	saveAs("Results", hsp_coloc_trajectories + Imagetitle + "_Hsp_colocal_traj.csv");
        run("Close");
selectWindow("Hsp.tif");
        run("Close");
roiManager("reset");

//filter Non-colocal tables
noncoloc_output_folder = trajectories_output + "/non-coloc/"
print(noncoloc_output_folder);
File.makeDirectory(noncoloc_output_folder);

//make subfolders for hsp vs client
client_noncoloc_trajectories = noncoloc_output_folder + "/Client/"
print(client_noncoloc_trajectories);
File.makeDirectory(client_noncoloc_trajectories);

hsp_noncoloc_trajectories = noncoloc_output_folder + "/Hsp/"
print(hsp_noncoloc_trajectories);
File.makeDirectory(hsp_noncoloc_trajectories);

//client
open(client_results);
IJ.renameResults("client");
open(Hsp_results);	
IJ.renameResults("Hsp");
run("count colocalized peaks Andrew", "table_1=client table_2=Hsp maximum_distance=3");
selectWindow("Hsp");
		run("Close")
selectWindow("client");
IJ.renameResults("Results");
runMacro(//INSERT PATH HERE:"/fiji-win32/Fiji.app/macros/Filter_for_noncolocal.ijm")
selectWindow("Colocalization");
close()
IJ.renameResults("client_non-colocalisation")
selectWindow("client_non-colocalisation")
saveAs("Results", output + "client_non-colocalisation.csv"); 
selectWindow("client_non-colocalisation.csv");
     run("Close"); 


open(client_image);
open("client_non-colocalisation.csv");
IJ.renameResults("Results")
runMacro(//INSERT PATH HERE:"/fiji-win32/Fiji.app/macros/results_to_ROI.txt"); // you need to put the location of your script
selectWindow("Results");
        run("Close");
runMacro(//INSERT PATH HERE:"/fiji-win32/Fiji.app/macros/peak_intensity.txt"); // you need to put the location of your script
selectWindow("Results");
        saveAs("Results", output + "client_non-colocal_traj.csv");
selectWindow("Results");
     	saveAs("Results", client_noncoloc_trajectories + Imagetitle + "_Client_non-colocal_traj.csv");
        run("Close");
selectWindow("client background corrected.tif");
        run("Close");
roiManager("reset");

//Hsp
open(Hsp_results);
IJ.renameResults("Hsp");
open(client_results);
IJ.renameResults("client");
run("count colocalized peaks Andrew", "table_1=Hsp table_2=client maximum_distance=3");

selectWindow("client");
		run("Close")
selectWindow("Hsp");
IJ.renameResults("Results");
runMacro(//INSERT PATH HERE:"/fiji-win32/Fiji.app/macros/Filter_for_noncolocal.ijm")
selectWindow("Colocalization");
		close()
IJ.renameResults("Hsp_non-colocalisation")
selectWindow("Hsp_non-colocalisation")
saveAs("Results", output + "Hsp_non-colocalisation.csv"); 
selectWindow("Hsp_non-colocalisation.csv");
     run("Close"); 

open(Hsp_image);	
open("Hsp_non-colocalisation.csv");
IJ.renameResults("Results")
runMacro(//INSERT PATH HERE:"/fiji-win32/Fiji.app/macros/results_to_ROI.txt"); // you need to put the location of your script
selectWindow("Results");
        run("Close");
runMacro(//INSERT PATH HERE:"/fiji-win32/Fiji.app/macros/peak_intensity.txt"); // you need to put the location of your script
//selectWindow("Results");
        saveAs("Results", output + "Hsp_non-colocal_traj.csv");
selectWindow("Results");
     	saveAs("Results", hsp_noncoloc_trajectories + Imagetitle + "_Hsp_non-colocal_traj.csv");
        run("Close");
selectWindow("Hsp.tif");
        run("Close");
roiManager("reset");
}
