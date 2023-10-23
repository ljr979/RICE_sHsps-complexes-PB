// compare all peaks (x, y positions) in table 1 with all peaks in table 2.

// when two peaks (not from the same table) are within a specified distance from each other

// only consider the two peaks that are closest to each other



// find peaks that are close to each other

// sort on distance

// pick top item (colocalized peaks) and remove these peaks from the other colocalized peaks

// put the output (distance) to table 1

// if the peaks are not colocalized leave the column empty

importClass(Packages.ij.IJ)
importClass(Packages.ij.WindowManager); 
importPackage(java.io); // import all java.io classes
   importClass(java.io.File) // import java.io.File
   importPackage(Packages.ij.text); // import ImageJ text package

importClass(Packages.ij.plugin.frame.RoiManager);
importClass(Packages.ij.gui.GenericDialog);

importPackage(Packages.ij.text);


 function print(s) {IJ.log(s);}

// default parameters

var maxDistance = 2.0;



// determine results tables

var nonImageWindows = WindowManager.getNonImageWindows();

var resultTables = [];

var titles = [];





for (var i = 0; i < nonImageWindows.length; i++) {

	

	if (nonImageWindows[i] instanceof TextWindow) {

	

		var textPanel = nonImageWindows[i].getTextPanel();

		var resultsTable = textPanel.getResultsTable();

		

		if (resultsTable) {

			resultTables.push(resultsTable);

			titles.push(nonImageWindows[i].getTitle());

		}

		

	}



}



if (titles.length < 2)

	IJ.error("This plugins needs at least 2 results tables.");



var dialog = new GenericDialog("Count Colocalized Peaks");

dialog.addChoice("table_1", titles, titles[0]);

dialog.addChoice("table_2", titles, titles[1]);

dialog.addNumericField("maximum_distance (in pixels)", maxDistance, 2);

dialog.showDialog();



if (dialog.wasCanceled())

	exit();

	

var index1 = dialog.getNextChoiceIndex();

var index2 = dialog.getNextChoiceIndex();

var table1 = resultTables[index1];

var table2 = resultTables[index2];

maxDistance = dialog.getNextNumber();



var colocalizations = new Array();



for (var i = 0; i < table1.getCounter(); i++) {



	var minDistance = maxDistance;

	var x1 = table1.getValue("X", i);

	var y1 = table1.getValue("Y", i);

	var s1 = table1.getValue("Frame", i);

	var colocalizedPeak = -1;

	

	for (var j = 0; j < table2.getCounter(); j++) {

	

		var x2 = table2.getValue("X", j);

		var y2 = table2.getValue("Y", j);

		var s2 = table2.getValue("Slice", j);

		

		var dx = x1 - x2;

		var dy = y1 - y2;

		var distance = Math.sqrt(dx * dx + dy * dy);

		

		if (distance < minDistance && s1 == s2) {

			minDistance = distance;

			colocalizedPeak = j;

		}

		

	}

	

	if (minDistance < maxDistance)

		colocalizations.push([i, colocalizedPeak, minDistance]);

}



// sort colocalization on distance

function compare(a, b) {

  return a[2] - b[2];

}



colocalizations.sort(compare);



for (var i = 0; i < table1.getCounter(); i++)

	table1.setValue("distance", i, -1);



while (colocalizations.length > 0) {



	var colocalization = colocalizations.pop();

	var row = colocalization[0];

	

	if (row != -1) {

		var x2 = table2.getValue("X", colocalization[1]);

		var y2 = table2.getValue("Y", colocalization[1]);

		var distance = colocalization[2];

		

		table1.setValue("X2", row, x2);

		table1.setValue("Y2", row, y2);

		table1.setValue("distance", row, distance);

		

		for (var i = colocalizations.length - 1; i >= 0; i--) {

			if (colocalizations[i][1] == colocalization[1])

				colocalizations[i][0] = -1;

		}

	}

}



table1.show(titles[index1]);



// draw colocalizations



// determine minimum and maximum

var xMin = table1.getValue("X", 0);

var xMax = table1.getValue("X", 0);

var yMin = table1.getValue("Y", 0);

var yMax = table1.getValue("Y", 0);

var sMin = table1.getValue("Frame", 0);

var sMax = table1.getValue("Frame", 0);



for (var i = 0; i < table1.getCounter(); i++) {

	var x = table1.getValue("X", i);

	var y = table1.getValue("Y", i);

	var s = table1.getValue("Frame", i);

	

	if (x < xMin) xMin = x;

	if (x > xMax) xMax = x;

	if (y < yMin) yMin = y;

	if (y > yMax) yMax = y;

	if (s < sMin) sMin = s;

	if (s > sMin) sMax = s;

}



for (var i = 0; i < table2.getCounter(); i++) {

	var x = table2.getValue("X", i);

	var y = table2.getValue("Y", i);

	var s = table2.getValue("Slice", i);

	

	if (x < xMin) xMin = x;

	if (x > xMax) xMax = x;

	if (y < yMin) yMin = y;

	if (y > yMax) yMax = y;

	if (s < sMin) sMin = s;

	if (s > sMin) sMax = s;

}



var slices = (sMax - sMin) + 1;



xMin = Math.floor(xMin);

xMax = Math.ceil(xMax);

yMin = Math.floor(yMin);

yMax = Math.ceil(yMax);



var width = xMax - xMin;

var height = yMax - yMin;



width *= 4;

height *= 4;



var image = IJ.createImage("Colocalization", "8-bit Black", width, height, slices);

var stack = image.getStack();



for (var i = 0; i < table1.getCounter(); i++) {

	var x = table1.getValue("X", i);

	var y = table1.getValue("Y", i);

	var s = (table1.getValue("Frame", i) - sMin) + 1;

	var distance = table1.getValue("distance", i);

	

	x -= xMin;

	y -= yMin;

	x *= 4;

	y *= 4;

	
	if (s <= stack.getSize()) {

		var ip = stack.getProcessor(s);

		ip.setColor(255);

		ip.drawOval(x - 2, y - 2, 4, 4);

	

		if (distance >= 0)

			ip.drawOval(x - maxDistance * 4, y - maxDistance * 4, maxDistance * 8, maxDistance * 8);
	}

}



for (var i = 0; i < table2.getCounter(); i++) {

	var x = table2.getValue("X", i);

	var y = table2.getValue("Y", i);

	var s = (table2.getValue("Slice", i) - sMin) + 1;

	

	x -= xMin;

	y -= yMin;

	x *= 4;

	y *= 4;

	
	if (s <= stack.getSize()) {

		var ip = stack.getProcessor(s);

		ip.setColor(255);

		ip.drawLine(x - 2, y - 2, x + 2, y + 2);

		ip.drawLine(x - 2, y + 2, x + 2, y - 2);
	}

}





image.getCalibration().xOrigin = -xMin * 4;

image.getCalibration().yOrigin = -yMin * 4;

image.getCalibration().pixelWidth = 0.25;

image.getCalibration().pixelHeight = 0.25;



image.show();
