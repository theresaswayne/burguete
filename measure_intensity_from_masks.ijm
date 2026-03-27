//@File(label = "Input directory", style = "directory") inputDir
//@File(label = "Output directory", style = "directory") outputDir
//@String (label = "File suffix", value = ".tif") fileSuffix
//@String (label = "Channel name to search for", value = "FUS") chanName
//@Integer(label = "Pixels to expand nucleus for defining cytoplasm", value = 3) dilate

// measure_intensity_from_masks.ijm
// ImageJ/Fiji script to measure a batch of images with corresponding masks, 1 cell per image

// *** Required inputs *** 
//  Masks and channel images must follow naming pattern:
//    (characters after hyphen are disregarded when finding the basename)
//  Images: XX_YYY-MaxIP_cell_ZZZ_DAPI.tif and XX_YYY-MaxIP_cell_ZZZ_FUS.tif
//  Masks: XX_YYY-MaxIP_cell_ZZZ_mask.tif and XX_YYY-MaxIP_cell_ZZZ_nucleus_mask.tif
//    where XX is sample, YYY is field, ZZZ is cell number

// Theresa Swayne, 2025
//  -------- Suggested text for acknowledgement -----------
//   "These studies used the Confocal and Specialized Microscopy Shared Resource 
//   of the Herbert Irving Comprehensive Cancer Center at Columbia University, 
//   funded in part through the NIH/NCI Cancer Center Support Grant P30CA013696."

// TO USE: Place all input images in the input folder. (not nested!). 
// 	Create a folder for the output files. 
//  Run the script in Fiji. 

// ---- Setup ----

close("*");

print("\\Clear"); // clear Log window

//setBatchMode(true); // faster performance but doesn't work for all functions
//run("Bio-Formats Macro Extensions"); // supports native microscope files

run("Set Measurements...", "area mean integrated stack display redirect=None decimal=3");
run("Roi Defaults...", "color=yellow stroke=2 group=0");

// collect data in a table with a time/date stamp
// collect data in a table with a time/date stamp
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
if (month<10) {monthString = "0"+(month+1);}
else {monthString = month+1;} // months start with 0
if (dayOfMonth<10) {dayString = +"0"+dayOfMonth;}
else {dayString = dayOfMonth;}
if (hour<10) {hourString = "0"+hour;}
else{hourString = hour;}
if (minute<10) {minString = "0" + minute;}
else{minString = minute;}
timeString = "" + year + "-" + monthString + "-" + dayString + "_" + hourString + "-" + minString; // have to start with empty string
dataName = timeString + "_Results.csv";
dataFile = outputDir + File.separator + dataName;
dataHeaders = "Filename,Cell_Area,Cell_Mean,Cell_IntDen,Cell_RawIntDen,Nucl_Area,Nucl_Mean,Nucl_IntDen,Nucl_RawIntDen,Cyto_Area,Cyto_Mean,Cyto_IntDen,Cyto_RawIntDen";
if (File.exists(dataFile)==false) { // start the file with headers
	File.append(dataHeaders, dataFile);	
	print("Created data file");
    }
// ---- Run ----

print("Starting");

processFolder(inputDir, outputDir, fileSuffix, chanName, dilate);

// Clean up images and get out of batch mode

close("*");
//setBatchMode(false);
setBatchMode("exit and display");
print("Finished");

// save Log
selectWindow("Log");
saveAs("text", outputDir + File.separator + timeString+ "_Log.txt");


// ---- Functions ----

function processFolder(input, output, suffix, chanName, dilate) {

	// this function searches for files matching the criteria and sends them to the processFile function
	filenum = -1;
	print("Processing folder", input);
	// scan folder tree to find files with correct suffix
	list = getFileList(input);
	list = Array.sort(list); // process in alphabetical order
	for (i = 0; i < list.length; i++) {
		//print("Checking file",i,"of",list.length);
		if(File.isDirectory(input + File.separator + list[i])) {
			processFolder(input + File.separator + list[i], output, suffix, chanName, dilate); // handle nested folders
		}
		if(startsWith(list[i], ".")) { // avoid dotfiles found on some mac disks
			continue;
		}
		if(endsWith(list[i], suffix)) {
			// check for channel name of interest
			candidateName = File.getNameWithoutExtension(list[i]);
			//print("Checking file", candidateName);
			if(candidateName.contains(chanName)) {
				filenum = filenum + 1;
				print("Found image file",filenum," at",candidateName);
				processFile(input, output, list[i], filenum, chanName, dilate); // pass the filename and parameters to the processFile function
			}
		}
	} // end of file list loop
} // end of processFolder function


function processFile(inputFolder, outputFolder, fileName, fileNumber, channelName, dilate) {
	
	// this function processes a single image
	
	// ---------- SETUP
	
	imagePath = inputFolder + File.separator + fileName;
	
	//print("Opening file",imagePath);
	// open the image file
	//run("Bio-Formats", "open=&imagePath");
	open(imagePath);
	// determine the name of the file without extension
	id = getImageID();
	dotIndex = lastIndexOf(fileName, ".");
	//basename = substring(fileName, 0, dotIndex);
	basename = File.getNameWithoutExtension(fileName);
	extension = substring(fileName, dotIndex);

	//print("Processing image",fileNumber," at path" ,imagePath);	
	
	selectImage(fileName);
	fluorID = getImageID(); // get the unique id of the image
	
	// open the corresponding masks
	lastUnderscore = lastIndexOf(basename, "_");
	cellName = substring(basename, 0, lastUnderscore);
	cellMaskName = cellName + "_mask.tif";
	cellMaskPath = inputFolder + File.separator + cellMaskName;
	nucMaskName = cellName + "_nucleus_mask.tif";
	nucMaskPath = inputFolder + File.separator + nucMaskName;

	if (!File.exists(cellMaskPath)) {
		print("Cell mask for",basename,"not found");
		close("*");
		return; // go to next file in folder
	}
	//print("Opening cell mask", cellMaskPath);
	open(cellMaskPath);
	selectWindow(cellMaskName);
	rename("cell");
	cellMaskID = getImageID();

	if (!File.exists(nucMaskPath)) {
		print("Nucleus mask for",basename,"not found");
		close("*");
		return; // go to next file in folder
	}
	//print("Opening nucleus mask", nucMaskPath);
	open(nucMaskPath);
	selectWindow(nucMaskName);
	rename("nucleus"); 
	nucMaskID = getImageID();
	
	// ---------- Measure and save
	
	run("Clear Results");
	
	// measure whole cell: first row of results
	
	selectImage("cell");
	while (!isActive(cellMaskID)) { // insist that the desired image is active
		wait(100);
	}
	run("Select None");
	getStatistics(area, mean, min, max);
	if(max == 0) {
		print("No segmented cell found for",fileName);
		close("*");
		return; // go to next file
	}
	setAutoThreshold("Default dark no-reset");
	//run("Threshold...");
	run("Create Selection");
	selectImage(fileName);
	while (!isActive(fluorID)) { // insist that the fluor image is active
		wait(100);
	}
	run("Restore Selection");
	run("Measure");
	run("Select None");
	//wait(100);
	
	// measure nucleus: 2nd row of results
	
	selectImage("nucleus");
	while (!isActive(nucMaskID)) { // insist that the desired image is active
		wait(100);
	}
	run("Select None");
	//wait(100);
	getStatistics(area, mean, min, max);
	if(max == 0) {
		print("No segmented nucleus found for",fileName);
		close("*");
		return; // go to next file
	}
	setAutoThreshold("Default dark no-reset");
	//run("Threshold...");
	run("Create Selection");
	selectImage(fileName);
	while (!isActive(fluorID)) { // insist that the fluor image is active
		wait(100);
	}
	run("Restore Selection");
	run("Measure");
	//wait(100);
	
	// dilate the nucleus and subtract its area from the cell
	
	selectImage("nucleus");
	while (!isActive(nucMaskID)) { // insist that the desired image is active
		wait(100);
	}
	run("Select None");
	setAutoThreshold("Default dark no-reset");
	//run("Threshold...");
	run("Create Selection");
	run("Enlarge...", "enlarge="+dilate+" pixel");
	//wait(100);
	selectImage("cell");
	while (!isActive(cellMaskID)) { // insist that the desired image is active
		wait(100);
	}
	run("Restore Selection");
	setBackgroundColor(0, 0, 0);
	run("Clear", "slice");
	run("Select None");
	setAutoThreshold("Default dark no-reset");
	//run("Threshold...");
	run("Create Selection");
	//wait(100);
	
	// measure the cytoplasm: 3rd row of results
	
	selectImage(fileName);
	while (!isActive(fluorID)) { // insist that the fluor image is active
		wait(100);
	}

	run("Restore Selection");
	run("Measure");
	run("Flatten");
	run("Scale...", "x=0.5 y=0.5 width=128 height=127 interpolation=Bilinear average create title=scaled");
	overlayName = basename + "_overlay.tif";
	saveAs("Tiff", outputFolder + File.separator + overlayName);
	close(overlayName);
	//wait(100);
	
	// Collect measurements 
	
	imageName = getResultString("Label", 0);
	resultsArray = newArray(imageName); 
	
	cellArea = getResult("Area", 0);
	cellMean = getResult("Mean", 0);
	cellIntDen = getResult("IntDen", 0);
	cellRawIntDen = getResult("RawIntDen", 0);
	
	resultsArray = Array.concat(resultsArray, cellArea, cellMean, cellIntDen, cellRawIntDen);
	
	nucArea = getResult("Area", 1);
	nucMean = getResult("Mean", 1);
	nucIntDen = getResult("IntDen", 1);
	nucRawIntDen = getResult("RawIntDen", 1);
	
	resultsArray = Array.concat(resultsArray, nucArea, nucMean, nucIntDen, nucRawIntDen);
	
	cytoArea = getResult("Area", 2);
	cytoMean = getResult("Mean", 2);
	cytoIntDen = getResult("IntDen", 2);
	cytoRawIntDen = getResult("RawIntDen", 2);
	
	resultsArray = Array.concat(resultsArray, cytoArea, cytoMean, cytoIntDen, cytoRawIntDen);
	
	resultsRow = String.join(resultsArray, ","); 
	
	print("Appending results for file",fileName);
	File.append(resultsRow, dataFile);	

	
	// ---------- CLEANUP
	
	close("*");
	//wait(1000);
	run("Clear Results");
	run("Collect Garbage");

} // end of processFile function


	