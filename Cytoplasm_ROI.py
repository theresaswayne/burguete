#@ File	(label = "Input directory", style = "directory") srcFile
#@ File	(label = "Output directory", style = "directory") dstFile
#@ String  (label = "File extension", value=".zip") ext
#@ String  (label = "File name contains", value = "RoiSet") containString
#@ Integer  (label = "ROIs to skip at beginning", value = "1") skipRois
#@ boolean (label = "Keep directory structure when saving", value = true) keepDirectories

# Cytoplasm_ROI.py
# Given an ROIset of a prescribed format, containing nucleis ad whole cell, 
# adds new ROIs for cytoplasm only
# ------- EXPECTED ROISet FORMAT ----------
# optional: 1st ROI = background measurement
# 2nd and 3rd: nucleus and whole cell, respectively, for cell # 1
# 4th and 5th, etc.: nucleus and cell for subsequent cells


from ij import IJ, ImagePlus, ImageStack
from ij.plugin import ZProjector
from ij.plugin.filter import RankFilters
# from ij.plugin.filter import BackgroundSubtracter
import net.imagej.ops
from net.imglib2.view import Views
from net.imglib2.img.display.imagej import ImageJFunctions as IL
from net.imglib2.algorithm.dog import DogDetection
from ij.gui import PointRoi
from jarray import zeros
from ij.measure import ResultsTable
from math import sqrt
from java.awt import Color
from ij.plugin.frame import RoiManager
from ij.gui import GenericDialog
import os
from loci.plugins import BF
from jarray import array

def process(srcDir, dstDir, currentDir, fileName, keepDirectories, skip):
	
 	IJ.run("Close All", "")
	
	imp = IJ.createImage("Dummy", "8-bit black", 2048, 2048, 1) # necessary to avoid RoiMgr errors
	imp.show()

	rm = RoiManager.getInstance()
	if not rm:
	  rm = RoiManager()
	rm.reset()
	
 	IJ.log("Processing ROI set:" + fileName)
 	
 	rm.runCommand("Open", os.path.join(currentDir, fileName))
 	
 	# loop through ROIs, skipping as needed
 	
 	numRois = rm.getCount()
 	IJ.log("This set has " + str(numRois) + " ROIs")
 	
 	# check for a valid number of ROIs
 	if (numRois - skipRois) % 2 != 0:
 		IJ.log("Skipped ROI set " + fileName + " because it has an invalid number of ROIs.")
 		return

	startIndex = skipRois # indices start at 0, so if we skip one, we start at 1
 	endIndex = numRois - startIndex # end is not included in range 
 	cellCount = 1
 	for RoiIndex in range(startIndex, endIndex, 2): # if we skip the first ROI, indices will be 1, 3, etc
		IJ.log("Processing ROI index " + str(RoiIndex))
 		rm.rename(RoiIndex, "Nucl_" + str(cellCount))
		rm.rename(RoiIndex+1, "Cell_" + str(cellCount))
		cellRois = array([RoiIndex, RoiIndex+1], 'i')
 		rm.setSelectedIndexes(cellRois)
 		IJ.log("Selected ROIs " + str(cellRois[0]) + "," + str(cellRois[1]))
 		rm.runCommand("XOR") # requires an image to be open
 		rm.runCommand("Add") # should be added at the end of the list
 		rm.deselect() # make sure nothing else selected
 		newTotal = rm.getCount()
 		IJ.log("There are now " + str(newTotal) + " ROIs")
 		rm.rename(newTotal-1, "Cyto_" + str(cellCount))
 		cellCount = cellCount + 1
 		
	# save the updated list
	saveDir = currentDir.replace(srcDir, dstDir) if keepDirectories else dstDir
	if not os.path.exists(saveDir):
		os.makedirs(saveDir)
	IJ.log("Saving to" + saveDir)
	rm.save(os.path.join(saveDir, fileName + "_Cyto_Rois.zip"))


def run():
	srcDir = srcFile.getAbsolutePath()
	dstDir = dstFile.getAbsolutePath()
	
	IJ.log("\\Clear")
	IJ.log("Processing ROIsets")
	
	for root, directories, filenames in os.walk(srcDir):
		filenames.sort();
	for filename in filenames:
		# Check for file extension
		if not filename.endswith(ext):
			continue
		# Check for file name pattern
		if containString not in filename:
			continue
		process(srcDir, dstDir, root, filename, keepDirectories, skipRois)

	rm = RoiManager.getInstance()
	if not rm:
		rm = RoiManager()
	rm.reset()
	IJ.run("Close All", "")
	IJ.log("Done!")

run()

#RoiManager.getName(index)