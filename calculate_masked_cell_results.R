# calculate_masked_cell_results.R

# Theresa Swayne, Columbia University, 2026
# -------- Suggested text for acknowledgement -----------
#   "These studies used the Confocal and Specialized Microscopy Shared Resource 
#   of the Herbert Irving Comprehensive Cancer Center at Columbia University, 
#   funded in part through the NIH/NCI Cancer Center Support Grant P30CA013696."

# --------- About this script ------------
# Takes a single table from the measure_intensity_from_masks.ijm script
#   and calculates cyto:nucleus and cyto:cell ratio for mean, IntDen, RawIntDen 
# Assumptions: The Filename column contains FUS image file names in the format:
#  XX_YYY-MaxIP_cell_ZZZ_FUS.tif
#    where XX is sample, YYY is field, ZZZ is cell number
# Output: CSV files, each row = 1 cell 
# Output data: original data plus:
#   Background-corrected mean, intden, rawintden for nucleus and cell
#   Cytoplasm area and Background-corrected mean/intden/rawintden
#   Cytoplasm/nucleus ratio (background corrected) for mean, intden

# ---- Setup and load data ----

require(tidyverse) # for data processing and plotting
require(stringr) # for string harvesting
require(tools) # for file name processing

# ---- Input and output setup ----

# Prompt for a file. No message will be displayed. Choose one of the files in the folder.
selectedFile <- file.choose()
inputFolder <- dirname(selectedFile) # the input is the parent of the selected file
fileName <- basename(selectedFile)

# Create an output folder with time-date stamp

#thisTime = format(Sys.time(),"%Y-%m-%d_%H%M")
#outputFolder <- file.path(inputFolder,paste0("Output_",thisTime))
#dir.create(outputFolder) # creates within the input folder if it does not already exist


# ---- Add sample info to the table ----

# parse image name
#  assuming format XX_YYY-MaxIP_cell_ZZZ_FUS.tif
#    where XX is sample, YYY is field, ZZZ is cell number
splitNames <- str_split(df$Filename, "_")

sampleNames <- sapply(splitNames, `[`, 1) # selects the 1st element of each vector in the list
fieldNames <- sapply(splitNames, `[`, 2)
fieldNames <- substr(fieldNames, 0,3) # remove -MaxIP
cellNames <- sapply(splitNames, `[`, 4)

df_mod <- df %>% mutate(
  `Sample` = sampleNames, 
  `Field` = fieldNames, 
  `Cell` = cellNames,
  .after = Filename
)

# ---- Calculate ratios ----

cytoNucMeanRatio <- df_mod$Cyto_Mean/df_mod$Nucl_Mean
cytoCellMeanRatio  <- df_mod$Cyto_Mean/df_mod$Cell_Mean

cytoNucIntDenRatio <- df_mod$Cyto_IntDen/df_mod$Nucl_IntDen
cytoCellIntDenRatio  <- df_mod$Cyto_IntDen/df_mod$Cell_IntDen

df_mod <- df_mod %>% mutate(
  `Mean_Cyto:Nucl` = cytoNucMeanRatio, 
  `Mean_Cyto:Cell` = cytoCellMeanRatio, 
  `IntDen_Cyto:Nucl` = cytoNucIntDenRatio, 
  `IntDen_Cyto:Cell` = cytoCellIntDenRatio, 
  .after = Cyto_RawIntDen
)

# ---- Summarize ----

df_summary <- df_mod %>% group_by(Sample) %>% summarise(
  `Cyto:Nuc Mean` = mean(`Mean_Cyto:Nucl`),
  `Cyto:Nuc IntDen` = mean(`IntDen_Cyto:Nucl`),
  `Cyto:Cell Mean` = mean(`Mean_Cyto:Cell`),
  `Cyto:Cell IntDen` = mean(`IntDen_Cyto:Cell`))

# ---- Visualize ----

p_cnm <-  ggplot(df_mod, aes(x=Sample,y=`Mean_Cyto:Nucl`, na.rm = TRUE)) +
  geom_boxplot() +
  scale_fill_manual(values = c("white", "grey")) +
  theme_minimal(base_size = 24) +
  stat_summary(fun=mean, geom="point", shape=18,
               size=3, color="red") +
  theme(legend.position="none",
        axis.title.y = element_text(size = 20)) +
  labs(y = "Cytoplasm:Nucleus Mean Intensity")

p_cni <-  ggplot(df_mod, aes(x=Sample,y=`IntDen_Cyto:Nucl`, na.rm = TRUE)) +
  geom_boxplot() +
  scale_fill_manual(values = c("white", "grey")) +
  theme_minimal(base_size = 24) +
  stat_summary(fun=mean, geom="point", shape=18,
               size=3, color="red") +
  theme(legend.position="none",
        axis.title.y = element_text(size = 20)) +
  labs(y = "Cytoplasm:Nucleus Integrated Intensity")

p_ccm <-  ggplot(df_mod, aes(x=Sample,y=`Mean_Cyto:Cell`, na.rm = TRUE)) +
  geom_boxplot() +
  scale_fill_manual(values = c("white", "grey")) +
  theme_minimal(base_size = 24) +
  stat_summary(fun=mean, geom="point", shape=18,
               size=3, color="red") +
  theme(legend.position="none",
        axis.title.y = element_text(size = 20)) +
  labs(y = "Cytoplasm:Cell Mean Intensity")

p_cci <-  ggplot(df_mod, aes(x=Sample,y=`IntDen_Cyto:Cell`, na.rm = TRUE)) +
  geom_boxplot() +
  scale_fill_manual(values = c("white", "grey")) +
  theme_minimal(base_size = 24) +
  stat_summary(fun=mean, geom="point", shape=18,
               size=3, color="red") +
  theme(legend.position="none",
        axis.title.y = element_text(size = 20)) +
  labs(y = "Cytoplasm:Cell Integrated Intensity")

# ---- Save results ----

outputName = paste0(fileName,"_calculations.csv")
write_csv(df_mod,file.path(inputFolder, outputName))

summaryName = paste0(fileName,"_summary.csv")
write_csv(df_summary,file.path(inputFolder, summaryName))

ggsave(file.path(inputFolder, paste0(fileName,"_cyto_nuc_mean_boxplot.png")), plot = p_cnm, width=7, height = 7)
ggsave(file.path(inputFolder, paste0(fileName,"_cyto_nuc_intden_boxplot.png")), plot = p_cni, width=7, height = 7)

ggsave(file.path(inputFolder, paste0(fileName,"_cyto_cell_mean_boxplot.png")), plot = p_ccm, width=7, height = 7)
ggsave(file.path(inputFolder, paste0(fileName,"_cyto_cell_intden_boxplot.png")), plot = p_cci, width=7, height = 7)

