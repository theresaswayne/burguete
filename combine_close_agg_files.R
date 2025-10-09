# combine_close_agg_files.R

# Theresa Swayne, Columbia University, 2025
# -------- Suggested text for acknowledgement -----------
#   "These studies used the Confocal and Specialized Microscopy Shared Resource 
#   of the Herbert Irving Comprehensive Cancer Center at Columbia University, 
#   funded in part through the NIH/NCI Cancer Center Support Grant P30CA013696."

# --------- About this script ------------
# Consolidates and performs calculations on a folder of ImageJ results tables in CSV format

# Important Assumptions for Input! 
#   Each spreadsheet is from a single 2D image containing one or more cells
#   Columns include Cell ID, Cytoplasm Area, Number of FUS Puncta in Cytoplasm
#   ImageJ is set to save Row numbers/Column numbers (Edit > Options > Input/Output)
# Output: One CSV file, each row = 1 cell 
# Output data: original data plus:
#   area in µm2 (based on pixel size in the pixSize variable)
#   puncta per unit area in each category

# ---- Setup and load data ----

require(tidyverse) # for data processing
require(stringr) # for string harvesting
require(tools) # for file name processing

# Parameters

# text to filter for in the end of the file name
finalText <- "Results.csv"
# width of a pixel in um
pixSize <- 0.1035718

# ---- Input and output setup ----

# Prompt for a file. No message will be displayed. Choose one of the files in the folder.
selectedFile <- file.choose()
inputFolder <- dirname(selectedFile) # the input is the parent of the selected file
outputFolder <- dirname(inputFolder) # parent of the input folder

# get file names
files <- dir(inputFolder, pattern = paste("*",finalText,sep=""))

# files <- list.files(inputFolder, pattern = "*.csv")


# tibble is used because of the warning that data_frame is deprecated.
mergedDataWithNames <- tibble(filename = files) %>% # tibble holding file names
  mutate(file_contents =
           map(filename,          # read files into a new data column
               ~ read_csv(file.path(inputFolder, .),
                          locale = locale(encoding = "latin1"),
                          na = c("", "N/A"))))

# sort by filename  
mergedDataWithNames <- mergedDataWithNames %>% 
  arrange(filename)

# unnest to make the list into a flat file again,
# but it now has 1 extra column to hold the filename
mergedDataFlat <- unnest(mergedDataWithNames, cols = c(file_contents))

# ---- Calculations ----

areaMicrons <- mergedDataFlat$`Cytoplasm Area` * (pixSize^2)
ch1PunctaPerSqUm <- mergedDataFlat$`FUS Puncta in Cytoplasm`/areaMicrons
ch2PunctaPerSqUm <- mergedDataFlat$`pFTAA Puncta in Cytoplasm`/areaMicrons

# both "close" columns should have identical numbers, or off by one. Pick the max
closeColumns <- select(mergedDataFlat, contains("within")) %>% 
  rowwise() %>% # group by rows to get the row-wise maximum
  mutate(Max = max(c_across(1:2))) # assumes there are only 2 columns

closePunctaPerSqUm <- closeColumns$Max/areaMicrons

# ---- Organize the data ----
  
# rename dummy "...1" column (ImageJ does not provide a header)
#data <- rename(data, "Measurement" = `...1`)
  
# sort by column 1 just to make sure we're in order (*may be unnecessary?)
#data <- arrange(data, Measurement)

outputData <- mergedDataFlat %>% 
  mutate("Cytoplasm Area um2" = areaMicrons, .after = "Cytoplasm Area") %>%
  mutate("FUS Puncta per um2" = ch1PunctaPerSqUm) %>%
  mutate("pFTAA Puncta per um2" = ch2PunctaPerSqUm) %>%
  mutate("Max Close Puncta per um2" = closePunctaPerSqUm)

# ---- Save results ----

outputFile = paste(basename(inputFolder), " merged", finalText, sep = "")
write_csv(outputData,file.path(outputFolder, outputFile))

