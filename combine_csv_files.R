# combine_csv_files.R
# merges a batch of csv files, putting the filename in a new column
# Requirements: All files must be within a single folder 
#	and the filenames to be processed must end with the pattern given in the finalText variable
# adapted from http://serialmentor.com/blog/2016/6/13/reading-and-combining-many-tidy-data-files-in-R
# Theresa Swayne, Columbia University, 2025
# -------- Suggested text for acknowledgement -----------
#   "These studies used the Confocal and Specialized Microscopy Shared Resource 
#   of the Herbert Irving Comprehensive Cancer Center at Columbia University, 
#   funded in part through the NIH/NCI Cancer Center Support Grant P30CA013696."

# Setup -------
require(tidyverse)
finalText <- "Results.csv"

# ---- Prompt for a data folder ----
# no message will be displayed. Choose one of the files in the folder
selectedFile <- file.choose()
inputFolder <- dirname(selectedFile) # the input is the parent of the selected file

# Read all the files in the folder ------

outputFolder <- dirname(inputFolder) # parent of the input folder

# get file names
files <- dir(inputFolder, pattern = paste("*",finalText,sep=""))

# tibble is used because of the warning that data_frame is deprecated.
mergedDataWithNames <- tibble(filename = files) %>% # tibble holding file names
  mutate(file_contents =
           map(filename,          # read files into a new data column
               ~ read_csv(file.path(inputFolder, .),
                          locale = locale(encoding = "latin1"),
                          na = c("", "N/A"))))

# unnest to make the list into a flat file again,
# but it now has 1 extra column to hold the filename
mergedDataFlat <- unnest(mergedDataWithNames, cols = c(file_contents))

# Write an output file of all the merged data ----------

outputFile = paste(basename(inputFolder), "merged", finalText) # spaces will be inserted
write_csv(mergedDataFlat,file.path(outputFolder, outputFile))

