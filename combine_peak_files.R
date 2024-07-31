# combine_peak_files.R
# imports, merges, sorts, and simplifies a batch of csv files exported from Volocity 6.3
# Requirement: All files must be within a single folder
# adapted from http://serialmentor.com/blog/2016/6/13/reading-and-combining-many-tidy-data-files-in-R

# Setup -------
require(tidyverse)

# ---- Prompt for a data folder ----
# no message will be displayed. Choose one of the files in the folder
selectedFile <- file.choose()
inputFolder <- dirname(selectedFile) # the input is the parent of the selected file

# Read all the files in the folder ------

outputFolder <- dirname(inputFolder) # parent of the input folder

# get file names
files <- dir(inputFolder, pattern = "*Results.csv") # avoid getting log files

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


# Correct the counts because 0,0 is always detected -------
# subtract 1 from columns 2-5
#mergedDataFlat[2-5] <- mergedDataFlat[2-5] - 3
mergedDataFlat<- mergedDataFlat %>% mutate(across(2:5,~.-1))

# Write an output file of all the merged data ----------

outputFile = paste(basename(inputFolder), "merged.csv") # spaces will be inserted
write_csv(mergedDataFlat,file.path(outputFolder, outputFile))

