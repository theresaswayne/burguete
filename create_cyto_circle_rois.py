#@ File	(label = "Input directory", style = "directory") srcFile
#@ File	(label = "Output directory", style = "directory") dstFile
#@ String  (label = "File extension", value=".zip") ext
#@ String  (label = "File name contains", value = "RoiSet") containString
#@ Boolean  (label = "Is there a background ROI at the beginning?", value = "True") BackgroundRoi
#@ Integer (label = "# Pixels to dilate nucleus before creating cytoplasm", value = "3") dilate
#@ Double (label = "Pixel size in microns", value = 0.1035718) pixSize
#@ Double (label = "Radius (um) for cytoplasmic circle", value = 13) radius
#@ Boolean (label = "Keep directory structure when saving", value = true) keepDirectories

# create_cyto_circle_rois.py
# Given an ROIset of a prescribed format, containing (optional: background) nucleus and whole cell, 
# adds new ROIs for cytoplasm only and for a cytoplasm limited by a desired radius
# and renames ROIs by compartment in cell number 

# ------- EXPECTED ROISet FORMAT ----------
# optional: 1st ROI = background measurement
# 2nd and 3rd: nucleus and whole cell, respectively, for cell # 1
# 4th and 5th, etc.: nucleus and cell for subsequent cells

# Theresa Swayne, 2025
#  -------- Suggested text for acknowledgement -----------
#   "These studies used the Confocal and Specialized Microscopy Shared Resource 
#   of the Herbert Irving Comprehensive Cancer Center at Columbia University, 
#   funded in part through the NIH/NCI Cancer Center Support Grant P30CA013696."

from ij import IJ, ImagePlus, ImageStack
from ij.plugin import ZProjector
from ij.plugin.filter import RankFilters
# from ij.plugin.filter import BackgroundSubtracter
import net.imagej.ops
from net.imglib2.view import Views
from net.imglib2.img.display.imagej import ImageJFunctions as IL
from ij.process import ImageStatistics as IS
from net.imglib2.algorithm.dog import DogDetection
from ij.gui import Roi, PointRoi, OvalRoi
from jarray import zeros
from ij.measure import ResultsTable
from ij.measure import Measurements
from ij.plugin.filter import Analyzer
from math import sqrt
from java.awt import Color
from ij.plugin.frame import RoiManager
from ij.gui import GenericDialog
import os
from loci.plugins import BF
from jarray import array
from ij.measure import ResultsTable

def process(imp, srcDir, dstDir, currentDir, fileName, keepDirectories, skip, dilate, pixSize, radius):
	
	#IJ.run("Close All", "")
	
	# Reference the current (dummy) image -- necessary to avoid RoiMgr errors
	imp = IJ.getImage() 

	# Create an ROI manager insance
	rm = RoiManager.getInstance()
	if not rm:
	  rm = RoiManager()
	rm.reset()
	
	IJ.log("Processing ROI set:" + fileName)
	
	# Open the ROI set
	rm.runCommand("Open", os.path.join(currentDir, fileName))
	
	# Check for a valid number of ROIs
	numRois = rm.getCount()

	IJ.log("This set has " + str(numRois) + " ROIs")
	# --- we expect 2 ROIs per cell, plus optional skipped ROIs. Skip the set if the count is wrong
	if (numRois - skip) % 2 != 0:
		IJ.log("Skipped ROI set " + fileName + " because it has an invalid number of ROIs.")
		return # exit the process file function

	# Create a cytoplasm ROI for each cell -- Loop through the cells via the ROI indices of the nuclear ROI
	
	startIndex = skip # indices start at 0, so if we skip one, we start at 1
	endIndex = numRois - startIndex # end is not included in range 
	cellIndex = 1 # this is the cell number we're working on
	for RoiIndex in range(startIndex, endIndex, 2): # if we skip the first ROI, indices will be 1, 3, etc
		IJ.log("Processing ROI index " + str(RoiIndex))
		
		# determine the center of the nucleus (before dilating) -- note these values are in scaled units
		rm.select(RoiIndex)
		stats = imp.getAllStatistics()
		nucX = stats.xCentroid
		nucY = stats.yCentroid
		IJ.log("The centroid of ROI " + str(RoiIndex) + " is " + str(nucX) + "," + str(nucY) + " um")
		# convert to pixels for circle construction
		cal = imp.getCalibration()
		pixelSize = cal.pixelWidth
		nucXpix = nucX/pixelSize
		nucYpix = nucY/pixelSize
		IJ.log("The centroid of ROI " + str(RoiIndex) + " is " + str(nucXpix) + "," + str(nucYpix) + " pixels")
		
		# dilate the nucleus ROI before defining cytoplasm, to avoid including bright nuclear fluorescence -- note that we use pixel units
		IJ.run("Enlarge...", "enlarge="+str(dilate) + " pixel")
		rm.runCommand(imp,"Update") 
		
		# rename the nucleus and cell ROIs with cellIndex
		rm.rename(RoiIndex, "Nucl_" + str(cellIndex)) #
		rm.rename(RoiIndex+1, "Cell_" + str(cellIndex))
		
		# define the cytoplasm as the cell exclusive of the nucleus 
		cellRois = array([RoiIndex, RoiIndex+1], 'i')
		rm.setSelectedIndexes(cellRois)
		IJ.log("Selected ROIs " + str(cellRois[0]) + "," + str(cellRois[1]))
		rm.runCommand("XOR") # requires an image to be open
		rm.runCommand("Add") # the new ROI should be added at the end of the list
		newTotal = rm.getCount()
		rm.deselect()
		rm.rename(newTotal-1, "Cyto_" + str(cellIndex))
		
		# define a cytoplasmic area constrained by a circle
		rm.deselect()
		# create a circle defined by the upper left corner, width, and height
		radiusPix = radius/pixelSize
		circleX = nucXpix-radiusPix
		circleY = nucYpix-radiusPix
		circle = OvalRoi(circleX, circleY, radiusPix*2, radiusPix*2)
		imp.setRoi(circle)
		rm.runCommand("Add") # circle is added at the end of the list
		newTotal = rm.getCount()
		rm.deselect()
		CytoCircleRois = array([newTotal-2, newTotal-1], "i") # select cyto and circle
		rm.setSelectedIndexes(CytoCircleRois)
		rm.runCommand("AND") # Create the intersection of the circle with the cytoplasm
		# Check to ensure that the intersection exists
		if imp.getRoi() == None:
			IJ.log("The circle and the cytoplasm do not intersect!")
		else:
			intStats = imp.getStatistics()
			IJ.log("The area of the intersection is " + str(intStats.area))
			rm.runCommand("Add") # the new ROI should be added at the end of the list
			rm.deselect()
			rm.rename(newTotal, "CytoConstrained_" + str(cellIndex)) # rename 
			# delete the circle
			rm.deselect()
			rm.select(newTotal-1) # select the circle
			rm.runCommand("Delete")
		
		# restore the original nucleus size for future measurement -- note that we use pixel units here
		rm.select(RoiIndex)
		IJ.run("Enlarge...", "enlarge="+str(-dilate)+" pixel")
		rm.runCommand(imp,"Update") 
		rm.deselect()
		
		cellIndex = cellIndex + 1
	
	numRois = rm.getCount()
	
	for RoiIndex in range(0, numRois):
		
		# add a line to the results table
		#table.incrementCounter()
		#table.addValue("Filename", fileName)
		
		rm.select(RoiIndex)
		roi = imp.getRoi()
		roiName = roi.getName()
		IJ.log("The name of ROI number " + str(RoiIndex) + " is " + roiName)
		# desiredName = "Cyto_"+roiNum;
		# if (matches(roiName, "Cyto_[0-9]{1,2}")) {
		#table.addValue("ROI name", roiName)
		#rm.runCommand(imp,"Measure");
		stats = imp.getStatistics(IS.AREA)
		#IJ.log("area: %s" %(stats.area))
		
		# Add to results table
		#table.addValue("ROI area",stats.area)
		rm.deselect() # make sure nothing else selected
		
	# rename the background ROI if present
	if skip == 1:
		rm.deselect()
		rm.rename(0, "Background")
		rm.deselect()
	
	# save the updated ROIs and area table
	saveDir = currentDir.replace(srcDir, dstDir) if keepDirectories else dstDir
	if not os.path.exists(saveDir):
		os.makedirs(saveDir)
	IJ.log("Saving ROIs to" + saveDir)
	baseName = os.path.splitext(fileName)[0]
	rm.deselect() # make sure nothing else selected
	rm.save(os.path.join(saveDir, baseName + "_Cyto_Rois.zip"))
	# table.save(os.path.join(saveDir, baseName + "_RoiAreas.csv"))

def run():
	srcDir = srcFile.getAbsolutePath()
	dstDir = dstFile.getAbsolutePath()
	
	IJ.log("\\Clear")
	IJ.run("Close All", "")
	IJ.run("Clear Results")
	IJ.log("Starting")
	
	imp = IJ.createImage("Dummy", "8-bit black", 2048, 2048, 1) # necessary to avoid RoiMgr errors
	imp.show()
	
	# define what to do if there is a background ROI
	if BackgroundRoi:
		skipRois = 1
	else:
		skipRois = 0
		
	# Set the scale for the dummy image to match the scale entered by the user
	
	origcal = imp.getCalibration()
	origPixSize = origcal.pixelWidth
	origUnit = origcal.getUnit()
	IJ.log("The original pixel size is " + str(origPixSize) + " " + origUnit)
	
	newcal = origcal # duplicate the existing calibration
	newcal.pixelWidth = pixSize # update the pixel size
	newcal.pixelHeight = pixSize
	newcal.setUnit("um")
	imp.setCalibration(newcal)
	
	# verify scale has been changed
	checkcal = imp.getCalibration()
	checkPixSize = checkcal.pixelWidth
	checkUnit = checkcal.getUnit()
	IJ.log("The new pixel size is " + str(checkPixSize) + " " + checkUnit)
	
	#table = ResultsTable()
	
	# Find files to process
	for root, directories, filenames in os.walk(srcDir):
		filenames.sort();
	for filename in filenames:
		IJ.log("Checking file " + filename)
		# Check for dotfile
		if filename.startswith("."):
			continue
		# Check for file extension
		if not filename.endswith(ext):
			continue
		# Check for file name pattern
		if containString not in filename:
			continue
		# Actually process the files
		process(imp, srcDir, dstDir, root, filename, keepDirectories, skipRois, dilate, pixSize, radius)

	rm = RoiManager.getInstance()
	if not rm:
		rm = RoiManager()
	rm.reset()
	
	#table.save(os.path.join(dstDir, "RoiAreas.csv"))
	
	IJ.run("Close All", "")
	IJ.run("Clear Results")
	IJ.log("Done!")

run()
