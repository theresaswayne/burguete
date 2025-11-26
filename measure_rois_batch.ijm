//@File(label = "Input image directory", style = "directory") inputDir
//@File(label = "Input ROIset directory", style = "directory") roiDir
//@File(label = "Output directory", style = "directory") outputDir
//@String (label = "File suffix", value = ".nd2") fileSuffix
//@Integer(label = "Channel to measure",value = 3) chan

// measure_rois_batch.ijm
// ImageJ/Fiji script to measure a batch of images with corresponding ROIsets
// One channel is measured (selected by user at runtime)

// *** Required inputs *** 
//  ROIset must be a Zip file with the same base name. Image is expected to end with -MaxIP 
//    (characters after hyphen are disregarded when finding the basename)
//  Image: No1_001-MaxIP
//  Roiset: no1_001_RoiSet_Cyto_Rois.zip

// Theresa Swayne, 2025
//  -------- Suggested text for acknowledgement -----------
//   "These studies used the Confocal and Specialized Microscopy Shared Resource 
//   of the Herbert Irving Comprehensive Cancer Center at Columbia University, 
//   funded in part through the NIH/NCI Cancer Center Support Grant P30CA013696."

// TO USE: Place all input images in the input folder. ROI sets can be in a different folder (not nested!). 
// 	Create a folder for the output files. 
// 	Place your desired processing steps in the processFile function.
// 	Collect any desired parameters in the script parameters at the top. 
//		See ImageJ wiki for more script parameter options.
//		Remember to pass your parameters into the processFolder and processFile functions!
//  Run the script in Fiji. 
//	Limitations -- cannot have >1 dots in the filename. 
//		Cannot have > 9 Cells per file.
// 		Each image must have an ROIset.

// ---- Setup ----

while (nImages>0) { // clean up open images
	selectImage(nImages);
	close();
}
print("\\Clear"); // clear Log window

roiManager("reset");

setBatchMode(true); // faster performance but doesn't work for all functions
run("Bio-Formats Macro Extensions"); // supports native microscope files

run("Set Measurements...", "area mean centroid integrated stack display redirect=None decimal=3");

// ---- Run ----

print("Starting");

processFolder(inputDir, roiDir, outputDir, fileSuffix, chan);

// Clean up images and get out of batch mode

while (nImages > 0) { // clean up open images
	selectImage(nImages);
	close(); 
}
//setBatchMode(false);
setBatchMode("exit and display");
print("Finished");


// ---- Functions ----

function processFolder(input, roiInput, output, suffix, channel) {

	// this function searches for files matching the criteria and sends them to the processFile function
	filenum = -1;
	print("Processing folder", input);
	// scan folder tree to find files with correct suffix
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i])) {
			processFolder(input + File.separator + list[i], output, suffix); // handle nested folders
		}
		if(startsWith(list[i], ".")) { // avoid dotfiles found on some mac disks
			continue;
		}
		if(endsWith(list[i], suffix)) {
		
			filenum = filenum + 1;
			processFile(input, roiInput, output, list[i], filenum); // pass the filename and parameters to the processFile function
		}
	} // end of file list loop
} // end of processFolder function


function processFile(inputFolder, roiFolder, outputFolder, fileName, fileNumber) {
	
	// this function processes a single image
	
	// ---------- SETUP
	
	imagePath = inputFolder + File.separator + fileName;
	
	// open the image file
	run("Bio-Formats", "open=&imagePath");
	
	// determine the name of the file without extension
	id = getImageID();
	dotIndex = lastIndexOf(fileName, ".");
	basename = substring(fileName, 0, dotIndex);
	extension = substring(fileName, dotIndex);

	print("Processing image",fileNumber," at path" ,imagePath, "with basename",basename);	
	
	// open the corresponding ROIset
	filenameParsed = split(basename, "-");
	roiFile = filenameParsed[0] + "_RoiSet_Cyto_Rois.zip"; // omit the -MaxIP and insert the ROI filename pattern
	roiPath = roiFolder + File.separator +roiFile;
	
	// open ROI
	roiManager("reset");
	print("Opening ROI", roiPath);
	roiManager("Open", roiPath);
	
	//numROIs = roiManager("count");	
	
	// ---------- Measure and save
	
	// make sure nothing is selected to begin with
	selectImage(id);
	roiManager("Deselect");
	
	// de-couple the rois from the channels in which they were drawn
	roiManager("Remove Channel Info");
	roiManager("Remove Slice Info");
	run("Select None");
	
//	for(roiIndex=0; roiIndex < numROIs; roiIndex++) { // loop through ROIs
//
//		selectImage(id);
//		roiNum = roiIndex + 1; // so that image names start with 1 like the ROI labels
//		roiManager("Select", roiIndex);  // ROI indices start with 0
//		roiName = Roi.getName();
//		//print("The name of ROI number",roiNum, "is",roiName);
//		//desiredName = "Cyto_"+roiNum;
//		if (matches(roiName, "Cyto_[0-9]")) {
//			// later: add circle generation
//		}
//	} // end of ROI loop

	// make sure nothing is selected to begin with
	selectImage(id);
	roiManager("Deselect");
	run("Select None");
	
	// do the measurement of all rois
	
	// option 1
	getDimensions(width, height, channels, slices, frames);
	if (channel > channels) {
		print("Channel selection is invalid for image",basename);
		return; // exit process file loop
	}
	Stack.setChannel(channel);
	print("Measuring channel",channel);
	roiManager("measure"); // measures all rois, one line per measurement, selected channel
	
	// option 2
	// roiManager("multi-measure measure_all"); // measures all channels -- first all rois on c1, then on c2, etc

	// option 3 -- must cycle through ROIs to make this work. Will be c1, c2, c3... for ROI1, etc.
	// run("Measure Stack...");
	
	// Save the table of measurements
	resultsName = basename + "_results.csv";
	selectWindow("Results");
	print("Saving results to",outputFolder+File.separator+resultsName);
	saveAs("Results", outputFolder+File.separator+resultsName);
	
	
	// ---------- CLEANUP
	
	run("Select None");
	//print("Saved",numROIs,"cropped ROIs.");
	selectImage(id);
	close();
	roiManager("Reset");
	run("Clear Results");
	run("Collect Garbage");


} // end of processFile function


	