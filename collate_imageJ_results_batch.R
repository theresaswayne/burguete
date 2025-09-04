# collate_imageJ_results_batch.R

# Theresa Swayne, Columbia University, 2025
# -------- Suggested text for acknowledgement -----------
#   "These studies used the Confocal and Specialized Microscopy Shared Resource 
#   of the Herbert Irving Comprehensive Cancer Center at Columbia University, 
#   funded in part through the NIH/NCI Cancer Center Support Grant P30CA013696."

# --------- About this script ------------
# calculates data for a folder of ImageJ results tables
# corrects for background based on an ROI measurement  
#   and calculates cytoplasm (cell - nucleus) mean, IntDen, RawIntDen 
#   from cells and nuclei measured in an ImageJ results table
# Important Assumptions for Input! 
#   Each spreadsheet is from a single 2D image containing one or more cells
#   At least Label, Mean, IntDen, RawIntDen, Area are included
#   ImageJ is set to save Row numbers/Column numbers (Edit > Options > Input/Output)
#   1st row is background measurement
#   subsequent rows are nuclei and cell, in order
#   Only the channel of interest is measured
# Output: Set of CSV files, each file = 1 image; each row = 1 cell 
# Output data: original data plus:
#   Background-corrected mean, intden, rawintden for nucleus and cell
#   Cytoplasm area and Background-corrected mean/intden/rawintden
#   Cytoplasm/nucleus ratio (background corrected) for mean, intden
#   Cytoplasm/cell ratio (background corrected) for mean, intden

# ---- Setup and load data ----

require(tidyverse) # for data processing
require(stringr) # for string harvesting
require(tools) # for file name processing

# ---- Input and output setup ----

# Prompt for a file. No message will be displayed. Choose one of the files in the folder.
selectedFile <- file.choose()
inputFolder <- dirname(selectedFile) # the input is the parent of the selected file

# Create an output folder with time-date stamp

thisTime = format(Sys.time(),"%Y-%m-%d_%H%M")
outputFolder <- file.path(inputFolder,paste0("Output_",thisTime))
dir.create(outputFolder) # creates within the input folder if it does not already exist

# Get names of CSV files in the folder

files <- list.files(inputFolder, pattern = "*.csv")

# ----- Function to process a single file ------

process_file_func <- function(f, out) {
  
  # read the file
  data <- read_csv(file.path(inputFolder, f), show_col_types = FALSE)
  
  # check for missing/extra data by counting rows
  if (nrow(data)%%2 == 0) {
    print(paste0(file, " has the wrong number of rows"))
    return()
  }
  
  # ---- Collect image info ----
  
  # get image name from first row in Label column, 
  #   from beginning to the first colon after a dot and 3 characters (file extension)
  measName <- data$Label[1]
  imageName <- str_match(measName, "^(.*)\\..{3}:")[2] # the second match is the one we want
  
  # get area per pixel, 
  #   by dividing the average IntDen by the average RawIntDen
  areaPerPixel <- mean(data$IntDen)/mean(data$RawIntDen)
  
  # get number of cells,
  #   assuming 2 measurements per cell and 1 row of background
  nCells <- (nrow(data) - 1)/2
  
  # get the background mean, 
  #   from the first data row
  background <- data$Mean[1]
  
  # ---- Organize the data ----
  
  # rename dummy "...1" column (ImageJ does not provide a header)
  data <- rename(data, "Measurement" = `...1`)
  
  # sort by column 1 just to make sure we're in order (*may be unnecessary?)
  data <- arrange(data, Measurement)
  
  # collect nuclei data
  #   collect the nuclei rows: 2, 4, 6, ...
  nuclei_rows <- seq(2, nrow(data), by=2)
  nuclei_data <- data[nuclei_rows,]
  
  #   add a cell number column to nuclei data
  nuclei_data <- mutate(nuclei_data, CellNumber=Measurement/2, .before = Measurement)
  
  #   rename all subsequent data columns
  nuclei_renamed <- rename_with(nuclei_data, 
                                ~ paste0("Nuc_", .x),
                                        any_of(c("Measurement", "Label", "Area", 
                                                 "Mean", "X", "Y", 
                                                 "IntDen", "RawIntDen")))
  
  # collect cell data
  #   collect the cell rows: 3, 5, 7...
  cells_rows <- seq(3, nrow(data), by=2)
  cells_data <- data[cells_rows,]
  
  #   add a cell number column to nuclei data
  cells_data <- mutate(cells_data, CellNumber=(Measurement-1)/2, .before = Measurement)
  
  #   rename all subsequent data columns
  cells_renamed <- rename_with(cells_data, 
                                ~ paste0("Cell_", .x),
                                any_of(c("Measurement", "Label", "Area", 
                                         "Mean", "X", "Y", 
                                         "IntDen", "RawIntDen")))
  
  
  # Merge the nuclei and cell data by the cell number, and sort by cell number
  data_merged <- full_join(nuclei_renamed, cells_renamed, by=join_by(CellNumber))
  
  # ---- Calculate results ----
  
  # Do background correction on nuclei and cells
  
  #   For Mean, subtract background value
  data_corrected <- data_merged %>% 
    mutate(Nuc_Corr_Mean = Nuc_Mean-background, .after=Nuc_Mean) %>%
    mutate(Cell_Corr_Mean = Cell_Mean-background, .after=Cell_Mean) 
  
  #   For IntDens, subtract background * area (um2)
  data_corrected <- data_corrected %>% 
    mutate(Nuc_Corr_IntDen = Nuc_IntDen-(background*Nuc_Area), .after=Nuc_IntDen) %>%
    mutate(Cell_Corr_IntDen = Cell_IntDen-(background*Cell_Area), .after=Cell_IntDen)
  
  #   For RawIntDens, subtract background * #pixels, where #pixels = area(um2)/area per pixel
  data_corrected <- data_corrected %>% 
    mutate(Nuc_Corr_RawIntDen = Nuc_RawIntDen-(background*Nuc_Area/areaPerPixel), .after=Nuc_RawIntDen) %>%
    mutate(Cell_Corr_RawIntDen = Cell_RawIntDen-(background*Cell_Area/areaPerPixel), .after=Cell_RawIntDen)
  
  # Calculate cytoplasm data
  
  data_calculated <- data_corrected %>%
    mutate(Cyto_Area = Cell_Area - Nuc_Area, .after=Cell_Corr_RawIntDen) %>%
    mutate(Cyto_Corr_RawIntDen = Cell_RawIntDen - Nuc_Corr_RawIntDen, .after=Cyto_Area)
  
  data_calculated <- data_calculated %>%
    mutate(Cyto_Corr_IntDen = Cyto_Corr_RawIntDen * areaPerPixel, .before = Cyto_Corr_RawIntDen) %>%
    mutate(Cyto_Corr_Mean = Cyto_Corr_IntDen/Cyto_Area, .after = Cyto_Area)
  
  # Calculate ratios
  #   cytoplasm/nucleus ratio (background corrected) for mean, intden, rawintden
  #   cytoplasm/cell ratio (background corrected) for mean, intden, rawintden
  # Corr_Cyto/Nuc_Mean,etc
  # Corr_Cyto/Cell_Mean, etc
  
  data_ratios <- data_calculated %>% 
    mutate(Ratio_CytoNuc_Corr_Mean = Cyto_Corr_Mean/Nuc_Corr_Mean) %>%
    mutate(Ratio_CytoNuc_Corr_IntDen = Cyto_Corr_IntDen/Nuc_Corr_IntDen) %>%
    mutate(Ratio_CytoCell_Corr_Mean = Cyto_Corr_Mean/Cell_Corr_Mean) %>%
    mutate(Ratio_CytoCell_Corr_IntDen = Cyto_Corr_IntDen/Cell_Corr_IntDen)
  
  # ---- Save results ----
  
  #outputFolder <- dirname(f) # parent of the input folder
  # generate filename from image name
  outputName = paste(imageName,"_calculations.csv", sep = "")
  # write CSV file
  write_csv(data_ratios,file.path(out, outputName))
  
  return()
} # end of process file function


# ---- Run the function on each file ----

for (file in files){
  process_file_func(file, outputFolder)
}

