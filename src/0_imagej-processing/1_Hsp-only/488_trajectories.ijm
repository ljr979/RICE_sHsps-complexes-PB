macro "AF488 CHANNEL background and peaks" 
 {

var path = File.openDialog("RAW image");
var dir = File.getParent(path);
var name = File.getName(path);

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

//ROI file is ROI.zip saved in the folder containing images. These are user defined ROIs that the macro calls on. after defining in IMAGEJ ROImanager, save as ROI.zip 
ROIpath = filesep + str + "ROI.zip";
print("Raw Image path= " + path)
print("ROI path= " + ROIpath);
output = dir + filesep
print("Output path= " + output);

run("Set Measurements...", "centroid stack redirect=None decimal=3");

blue_results = output + "blue_results.csv" // this is the path for where the table results are

blue_image = output + "blue background corrected.tif"


//saved as same ROI size, gaussian blur image, saved in folder same place as ROI.zip.
blue_background_image = filesep + str + "blue_background.tif"

//folder where trajectories get saved
trajectories_path = filesep + str + "/Trajectories/";
print(trajectories_path);

File.makeDirectory(trajectories_path);
open(path);	//open original raw image
Imagetitle = getTitle();
rename("Raw");

roiManager("Open", ROIpath);
roiManager("Select", 0);	//this will be the channel with molecules in it (named 0)

run("Duplicate...", "title=blue duplicate");
	
//Background correction - subtract gaussian blur, beam inequality
open(blue_background_image);
rename("AVG_blue");
selectWindow("blue");
run("32-bit");
run("Beam Profile Correction", "electronic_offset=1300 background_image=AVG_blue");
setMinAndMax(0, 65536);
run("16-bit");
saveAs("tiff", output + "blue background corrected.tif")


open(blue_image);
//crop second time in order to remove black edges in Hsp channel (this is slighty smaller than original ROI)
roiManager("Select", 1);
run("Duplicate...", "title=[blue background corrected] duplicate");
saveAs("tiff", output + "blue background corrected.tif");	//this will save the movie that can be used to get the intensity trajectories
//saves max projection of first 20 frames to use for identifying locations of molecules
run("Z Project...", "stop=20 projection=[Max Intensity]");
saveAs("tiff", output + "MAX_blue.tif");
roiManager("reset");
selectWindow("MAX_blue.tif");
run("Enhance Contrast", "saturated=0.35");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");


//Find the peaks within the image
run("Peak Finder");// need to check the threshold here
waitForUser("Have you selected your peaks? \n \nPress OK to continue");
roiManager("Measure");	//this will allow for a table to be produced
saveAs("Results", output + "blue_results.csv"); 
close();
selectWindow("Results"); 
     run("Close"); 
roiManager("reset");

open(blue_image);
open("blue_results.csv");
IJ.renameResults("Results")
//saves the locations of the peaks to ROIs so we can extract intensity
//runMacro("results_to_ROI.txt"); you need to put the location of your script (locally)
runMacro( "/Users/ljr979/Desktop/Fiji.app/plugins/Macros/results_to_ROI.txt");
selectWindow("Results");
        run("Close");
//extracts peak intensity at all these ROIs
runMacro("/Users/ljr979/Desktop/Fiji.app/plugins/Macros/peak_intensity.txt"); // you need to put the location of your script
selectWindow("Results");
//saves trajectory
        saveAs("Results", output + "blue_traj.csv");
        
selectWindow("Results");
        saveAs("Results",trajectories_path + Imagetitle + ".csv");
        run("Close");
selectWindow("blue background corrected.tif");
        run("Close");
selectWindow("blue background corrected-1.tif");
        run("Close");
selectWindow("blue background corrected-2.tif");
        run("Close");
roiManager("reset");
selectWindow("AVG_blue")
		run("Close")
selectWindow("Raw")
		run("Close")
roiManager("reset");


}
