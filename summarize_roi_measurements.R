# script (or commands) to calculate summaries for ImageJ data

# Theresa Swayne, Columbia University, 2024-25
# -------- Suggested text for acknowledgement -----------
#   "These studies used the Confocal and Specialized Microscopy Shared Resource 
#   of the Herbert Irving Comprehensive Cancer Center at Columbia University, 
#   funded in part through the NIH/NCI Cancer Center Support Grant P30CA013696."

# ---- Setup -------
require(tidyverse)
require(skimr)
require(ggplot2)

# ---- Load data ----
# no message will be displayed.
selectedFile <- file.choose()

df <- read_csv(selectedFile, show_col_types = FALSE)

# ---- Get an overview of the data ----
skim(df) # notes missing data
summary(df) # mean, median, min, max, quartiles
summary(df$Ratio_CytoConstrToCell_Corr_Mean) # for a specific column

# ---- Plot important variables ----
hist(df$Ratio_CytoConstrToNuc_Corr_Mean, xlim = c(0,1), labels = FALSE)
hist(df$Ratio_CytoToNuc_Corr_Mean, xlim = c(0,1))
#hist(data$Ratio_CytoNuc_Corr_Mean, xlim = c(0,1))
p <- ggplot(data = df, aes(x = Ratio_CytoConstrToNuc_Corr_Mean)) +
  geom_histogram() +
  xlim = c(0,1)
ggsave("foo.png")
p