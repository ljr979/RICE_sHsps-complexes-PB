macro "647 CHANNEL background and peaks" 
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
ROIpath = filesep + str + "ROI.zip";
print("Raw Image path= " + path)
print("ROI path= " + ROIpath);
output = dir + filesep
print("Output path= " + output);

run("Set Measurements...", "centroid stack redirect=None decimal=3");

client_results = output + "client_results.csv" // this is the path for where the table results are
client_image = output + "client background corrected.tif"
client_background_image = filesep + str + "client_background.tif" 

//creates folder to save trajectories in
trajectories_path = filesep + str + "/Trajectories/";

print(trajectories_path);

File.makeDirectory(trajectories_path);

open(path);	//open original raw image
Imagetitle = getTitle();
rename("Raw");
roiManager("Open", ROIpath);
roiManager("Select", 0);	//this will be the channel with molecules in it
run("Duplicate...", "title=client duplicate");
	

//Background correction (fixes any uneven aspects of the beam profile)(gaussian blur)
open(client_background_image);
rename("AVG_client");
selectWindow("client");
run("32-bit");
run("Beam Profile Correction", "electronic_offset=400 background_image=AVG_client");
setMinAndMax(0, 65536);
run("16-bit");
saveAs("tiff", output + "client background corrected.tif")



open(client_image);
//crop second time in order to remove black edges in channel (smaller than the previous one)
//define this ROI at same time as the other one.
roiManager("Select", 1);
run("Duplicate...", "title=[client background corrected] duplicate");
saveAs("tiff", output + "client background corrected.tif");	//this will save the movie that can be used to get the intensity trajectories
run("Z Project...", "stop=20 projection=[Max Intensity]");
saveAs("tiff", output + "MAX_client.tif");
roiManager("reset");
selectWindow("MAX_client.tif");
run("Enhance Contrast", "saturated=0.35");
run("Properties...", "channels=1 slices=1 frames=1 unit=pixels pixel_width=1 pixel_height=1 voxel_depth=1.0000000");

//Find the peaks within the image
run("Peak Finder");// need to check the threshold here
waitForUser("Have you selected your peaks? \n \nPress OK to continue");
roiManager("Measure");	//this will allow for a table to be produced
saveAs("Results", output + "client_results.csv"); 
close();
selectWindow("Results"); 
     run("Close"); 
roiManager("reset");

open(client_image);
open("client_results.csv");
IJ.renameResults("Results")
//runMacro("results_to_ROI.txt"); // you need to put the location of your script (local location of script)
runMacro( "/Users/ljr979/Desktop/Fiji.app/plugins/Macros/results_to_ROI.txt");
selectWindow("Results");
        run("Close");
        //extract trajectories
runMacro("/Users/ljr979/Desktop/Fiji.app/plugins/Macros/peak_intensity.txt"); // you need to put the location of your script
selectWindow("Results");
        saveAs("Results", output + "client_traj.csv");
        
selectWindow("Results");
        saveAs("Results",trajectories_path + Imagetitle + ".csv");
        run("Close");
selectWindow("client background corrected.tif");
        run("Close");
roiManager("reset");
selectWindow("AVG_client")
		run("Close")
selectWindow("Raw")
		run("Close")
roiManager("reset");


}
