# calculate_ImageJ_results_batch.R

# Theresa Swayne, Columbia University, 2025
# -------- Suggested text for acknowledgement -----------
#   "These studies used the Confocal and Specialized Microscopy Shared Resource 
#   of the Herbert Irving Comprehensive Cancer Center at Columbia University, 
#   funded in part through the NIH/NCI Cancer Center Support Grant P30CA013696."

# --------- About this script ------------
# calculates data for a folder of ImageJ results tables
# corrects for background based on an ROI measurement  
#   and calculates background corrected mean, IntDen, RawIntDen 
#   from compartments measured in an ImageJ results table
# Important Assumptions for Input! 
#   Each spreadsheet is from a single 2D image containing one or more cells
#   At least Label, Mean, IntDen, RawIntDen, Area are included
#   ImageJ is set to save Row numbers/Column numbers (Edit > Options > Input/Output)
#   1st row is background measurement
#   subsequent rows are nuclei, cell in order, followed by cytoplasm and constrained cytoplasm, each numbered by cell
#   Only the channel of interest is measured
# Output: Set of CSV files, each file = 1 image; each row = 1 cell 
# Output data: original data plus:
#   Background-corrected mean, intden, rawintden for nucleus and cell
#   Cytoplasm area and Background-corrected mean/intden/rawintden
#   Cytoplasm/nucleus ratio (background corrected) for mean, intden
#   Cytoplasm/cell ratio (background corrected) for mean, intden
# cytoplasm is calculated as the fll cytoplasm plus a constrained cytoplasm

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
# arguments are an input folder and file and an output folder

process_file_func <- function(input, f, out) {
  
  # read the file
  data <- read_csv(file.path(input, f), show_col_types = FALSE)
  
  # check for missing/extra data by counting rows
  # we should have 4 compartments per cell, plus one background per sheet 
  if ((nrow(data) - 1)%%4 != 0) {
    print(paste0(file, " has the wrong number of rows"))
    return()
  }
  
  # ---- Collect image info ----
  
  # get image name from first row in Label column, 
  #   from beginning to the first colon after a dot and 3 characters (file extension)
  #measName <- data$Label[1]
  #imageName <- str_match(measName, "^(.*)\\..{3}:")[2] # the second match is the one we want
  
  # get area per pixel, 
  #   by dividing the average IntDen by the average RawIntDen
  areaPerPixel <- mean(data$IntDen)/mean(data$RawIntDen)
  
  # get number of cells,
  #   assuming 4 measurements per cell and 1 row of background
  nCells <- (nrow(data) - 1)/4
  
  # get the background mean, 
  #   from the first data row
  background <- data$Mean[1]
  
  # ---- Clean the data ----
  
  # rename dummy "...1" column (ImageJ does not provide a header)
  data <- rename(data, "Measurement" = `...1`)
  
  # sort by column 1 just to make sure we're in order (*may be unnecessary?)
  data <- arrange(data, Measurement)
  
  # since the compartments are not stored in cell order, 
  #   we need to parse the Label to find the compartment and cell number
  data_cleaned <- 
    data %>% 
    separate_wider_delim(Label,
                         ":",
                         names = c("Image", "Roi", "c", "channelEtc"))
  
  # parse the channelEtc column to get the channel -- unnecessary if channel is included in results
  # data_parsed <- data_parsed %>% mutate(channel = substring(channelEtc, 1,1), .after=Roi)
  
  # parse the Roi column to get cell number (background will be NA)
  data_cleaned <- data_cleaned %>% separate_wider_delim(Roi, "_", names = c("Roi", "CellNumber"), too_few = "align_start")
  
  # remove the useless c and channelEtc columns
  data_cleaned <- data_cleaned %>% select(-c("c", "channelEtc"))
  
  # ---- Re-organize the data so each row represents one cell ----
  
  # collect nuclei data
  #   collect the nuclei rows: 2, 6, 10...
  #nuclei_rows <- seq(2, nrow(data), by=4)
  #nuclei_data <- data[nuclei_rows,]
  
  #   better: collect by ROI name
  nuclei_data <- filter(data_cleaned, Roi == "Nucl")

  #   rename all subsequent data columns
  nuclei_renamed <- rename_with(nuclei_data, 
                                ~ paste0("Nuc_", .x),
                                        any_of(c("Measurement", "Label", "Area", 
                                                 "Mean", "X", "Y", 
                                                 "IntDen", "RawIntDen", "Ch", "Roi")))
  
  # collect cell data
  cells_data <- filter(data_cleaned, Roi == "Cell")
  
  #   rename all subsequent data columns
  cells_renamed <- rename_with(cells_data, 
                                ~ paste0("Cell_", .x),
                                any_of(c("Measurement", "Label", "Area", 
                                         "Mean", "X", "Y", 
                                         "IntDen", "RawIntDen", "Ch", "Roi")))
  
  
  # collect cyto data
  cyto_data <- filter(data_cleaned, Roi == "Cyto")
  
  #   rename all subsequent data columns
  cyto_renamed <- rename_with(cyto_data, 
                               ~ paste0("Cyto_", .x),
                               any_of(c("Measurement", "Label", "Area", 
                                        "Mean", "X", "Y", 
                                        "IntDen", "RawIntDen", "Ch", "Roi")))
  
  # collect constrained cytoplasm data
  cytoConstr_data <- filter(data_cleaned, Roi == "CytoConstrained")
  
  #   rename all subsequent data columns
  cytoConstr_renamed <- rename_with(cytoConstr_data, 
                               ~ paste0("CytoConstr_", .x),
                               any_of(c("Measurement", "Label", "Area", 
                                        "Mean", "X", "Y", 
                                        "IntDen", "RawIntDen", "Ch", "Roi")))
  
  # Merge the data by the cell number, and don't duplicate the image name column
  data_merged <- full_join(nuclei_renamed, cells_renamed, by=join_by(CellNumber, Image))
  data_merged <- right_join(data_merged, cyto_renamed, by=join_by(CellNumber, Image))
  data_merged <- right_join(data_merged, cytoConstr_renamed, by=join_by(CellNumber, Image))
  
  # Add a background column
  data_merged <- data_merged %>% mutate(BackgroundMean = background, .after=CellNumber)
  
  
  # ---- Calculate results ----
  
  # Do background correction on each compartment
  
  #   For Mean, subtract background mean value
  data_corrected <- data_merged %>% 
    mutate(Nuc_Corr_Mean = Nuc_Mean-BackgroundMean, .after=Nuc_Mean) %>%
    mutate(Cell_Corr_Mean = Cell_Mean-BackgroundMean, .after=Cell_Mean) %>%
    mutate(Cyto_Corr_Mean = Cyto_Mean-BackgroundMean, .after=Cyto_Mean) %>%
    mutate(CytoConstr_Corr_Mean = CytoConstr_Mean-BackgroundMean, .after=CytoConstr_Mean)
  
  #   For IntDens, subtract background * area (um2)
  data_corrected <- data_corrected %>% 
    mutate(Nuc_Corr_IntDen = Nuc_IntDen-(BackgroundMean*Nuc_Area), .after=Nuc_IntDen) %>%
    mutate(Cell_Corr_IntDen = Cell_IntDen-(BackgroundMean*Cell_Area), .after=Cell_IntDen) %>%
    mutate(Cyto_Corr_IntDen = Cyto_IntDen-(BackgroundMean*Cyto_Area), .after=Cyto_IntDen) %>%
    mutate(CytoConstr_Corr_IntDen = CytoConstr_IntDen-(BackgroundMean*CytoConstr_Area), .after=CytoConstr_IntDen)
  
  #   For RawIntDens, subtract background * #pixels, where #pixels = area(um2)/area per pixel
  data_corrected <- data_corrected %>% 
    mutate(Nuc_Corr_RawIntDen = Nuc_RawIntDen-(BackgroundMean*Nuc_Area/areaPerPixel), .after=Nuc_RawIntDen) %>%
    mutate(Cell_Corr_RawIntDen = Cell_RawIntDen-(BackgroundMean*Cell_Area/areaPerPixel), .after=Cell_RawIntDen) %>%
    mutate(Cyto_Corr_RawIntDen = Cyto_RawIntDen-(BackgroundMean*Cyto_Area/areaPerPixel), .after=Cyto_RawIntDen) %>%
    mutate(CytoConstr_Corr_RawIntDen = CytoConstr_RawIntDen-(BackgroundMean*CytoConstr_Area/areaPerPixel), .after=CytoConstr_RawIntDen)

  # Calculate ratios
  #   cytoplasm/nucleus ratio (background corrected) for mean, intden, rawintden
  #   cytoplasm/cell ratio (background corrected) for mean, intden, rawintden
  # Corr_Cyto/Nuc_Mean,etc
  # Corr_Cyto/Cell_Mean, etc
  
  data_ratios <- data_calculated %>% 
    mutate(Ratio_CytoToNuc_Corr_Mean = Cyto_Corr_Mean/Nuc_Corr_Mean) %>%
    mutate(Ratio_CytoToNuc_Corr_IntDen = Cyto_Corr_IntDen/Nuc_Corr_IntDen) %>%
    mutate(Ratio_CytoToCell_Corr_Mean = Cyto_Corr_Mean/Cell_Corr_Mean) %>%
    mutate(Ratio_CytoToCell_Corr_IntDen = Cyto_Corr_IntDen/Cell_Corr_IntDen) %>%
     mutate(Ratio_CytoConstrToNuc_Corr_Mean = CytoConstr_Corr_Mean/Nuc_Corr_Mean) %>%
     mutate(Ratio_CytoConstrToNuc_Corr_IntDen = CytoConstr_Corr_IntDen/Nuc_Corr_IntDen) %>%
     mutate(Ratio_CytoConstrToCell_Corr_Mean = CytoConstr_Corr_Mean/Cell_Corr_Mean) %>%
     mutate(Ratio_CytoConstrToCell_Corr_IntDen = CytoConstr_Corr_IntDen/Cell_Corr_IntDen)
              
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
  process_file_func(inputFolder, file, outputFolder)
}

