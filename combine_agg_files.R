# combine_agg_files.R

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
#   puncta per unit area

# ---- Setup and load data ----

require(tidyverse) # for data processing
require(stringr) # for string harvesting
require(tools) # for file name processing

# Parameters

# text to filter for in the end of the file name
finalText <- "Results.csv"
# width of a pixel in um
pixSize <- 0.1

# ---- Input and output setup ----

# Prompt for a file. No message will be displayed. Choose one of the files in the folder.
selectedFile <- file.choose()
inputFolder <- dirname(selectedFile) # the input is the parent of the selected file

# Create an output folder with time-date stamp

thisTime = format(Sys.time(),"%Y-%m-%d_%H%M")
outputFolder <- file.path(inputFolder,paste0("Output_",thisTime))
dir.create(outputFolder) # creates within the input folder if it does not already exist

# Get names of CSV files in the folder


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

# sort by filename and add a column for the row number (timepoint) 
mergedDataWithNames <- mergedDataWithNames %>% 
  arrange(filename) %>%
  mutate(row_number = row_number())

# unnest to make the list into a flat file again,
# but it now has 1 extra column to hold the filename
mergedDataFlat <- unnest(mergedDataWithNames, cols = c(file_contents))

# ---- Calculations ----

areaMicrons <- mergedDataFlat$`Cytoplasm Area` * (pixSize^2)
punctaPerSqUm <- mergedDataFlat$`Number of FUS Puncta in Cytoplasm`/areaMicrons

# ---- Organize the data ----
  
# rename dummy "...1" column (ImageJ does not provide a header)
#data <- rename(data, "Measurement" = `...1`)
  
# sort by column 1 just to make sure we're in order (*may be unnecessary?)
#data <- arrange(data, Measurement)

outputData <- mergedDataFlat %>% 
  mutate("Cytoplasm Area um2" = areaMicrons, .after = "Cytoplasm Area") %>%
  mutate("Puncta per um2" = punctaPerSqUm, .after = "Number of FUS Puncta in Cytoplasm")

# ---- Save results ----

outputFile = paste(basename(inputFolder), " merged", finalText, sep = "")
write_csv(outputData,file.path(outputFolder, outputFile))




