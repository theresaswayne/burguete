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
summary(No6_1118$Ratio_CytoNuc_Corr_Mean) # mean, median, min, max, quartiles
summary(df$Ratio_CytoConstrToCell_Corr_Mean) # for a specific column

# ---- Plot important variables ----

bins <- seq (0,1,by=0.05)

# for new data
hist(No36_1202$Ratio_CytoConstrToNuc_Corr_Mean, xlim = c(0,1), breaks = bins, labels = FALSE)
hist(No36_1202$Ratio_CytoToNuc_Corr_Mean, breaks = bins, xlim = c(0,1))

# for old data
hist(No36_1118$Ratio_CytoNuc_Corr_Mean, breaks = 20, xlim = c(0,1))

p <- ggplot(data = No6, aes(x = Ratio_CytoConstrToNuc_Corr_Mean)) +
  geom_histogram() +
  xlim = c(0,1)

plt <- ggplot(data = diamonds) +
  geom_bar(mapping = aes(x = cut, y = after_stat(prop), group = 1))

plt <- ggplot(data = No6) +
  geom_bar(mapping = aes(x = Ratio_CytoConstrToNuc_Corr_Mean, 
                         y = after_stat(prop))) +

#ggsave("foo.png")
p