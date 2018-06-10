################# Identify Regions under Selection from Sliding Window Analysis #########
#
# This script will filter the output file from sliding window analysis...
#    to include only regions where the FST average goes above the confidence interval
# You can then identify which loci are within the selected regions
#
# MF 3/22/2018 for PCod Compare Project
#
#########################################################################################


# Install Packages --------------------------------------------------------
install.packages("dplyr")
install.packages("ggplot2")
library(readr)
library(ggplot2)
library(dplyr)



# Load Data -------------------------------------------------------------
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/SlidingWindow")
source("SlidingWindow_FindMarkers_Function.R")
source("Plot_SlidingWindowAnalysis_Functions.R")

mydata <- read_delim("EastvWest/batch_8_SWA_eastwest_globalFST_filtered_output_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150_FILTERED.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
head(mydata)




# Add columns for selection -------------------------------------------------------------
## OPTION 1: output only divergent regions
mydata_selection <- mydata %>%
  mutate(positive = ifelse(`Fst.Fct` > upper_95, "1", "0")) %>%
  mutate(negative = ifelse(`Fst.Fct` < lower_95, "1", "0")) %>%
  mutate(selection = ifelse(positive == 1, "Yes", "No")) 
head(mydata_selection)

## OPTION 2: output both divergent and balancing regions
mydata_selection <- mydata %>%
  mutate(positive = ifelse(`Fst.Fct` > upper_95, "1", "0")) %>%
  mutate(negative = ifelse(`Fst.Fct` < lower_95, "1", "0")) %>%
  mutate(selection = ifelse(positive == 1 | negative == 1, "Yes", "No")) 
head(mydata_selection)




# Visualize Fst v. Bootstrap CI ----------------------------------------------
## positive selection: Fst above upper 95% confidence interval
ggplot(data=mydata_selection, aes(x = `Fst.Fct`, y = upper_95)) +
  geom_point(aes(color=positive)) +
  labs(title="mydata Sliding Window\nPositive Selection")

## negative selection: Fst below lower 95% confidence interval
ggplot(data=mydata_selection, aes(x = `Fst.Fct`, y = lower_95)) +
  geom_point(aes(color=negative)) +
  labs(title="mydata Sliding Window\nNegative Selection")






# Visualize P values ------------------------------------------------------
ggplot(data=mydata_selection, aes(x = `Fst.Fct`, y = pvalue)) +
  geom_point(aes(color=selection)) +
  ylim(-0.1,1) +
  xlim(-0.1,1) +
  labs(title="mydata Sliding Window\nP values")






# Match Data sets on Plot ---------------------------------------------------------
## plot all chromosomes using overlay function
colnames(mydata_selection) <- c("chromosome","position","Fst/Fct","Mean_boostrap","lower_95","upper_95","pvalue","positive","negative","selection")



just_plot_diverge(data1 = mydata_selection, Nb_divisions = 150, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="plots/Eastwest/SWA_EastWest_Divergence_filtered")







# Filter Data and Write to File --------------------------------------------
## this pipeline keeps row numbers in the "window" column so that you can see which windows are next to each other. 
mydata_filter <- mydata_selection %>%
  tibble::rownames_to_column("Window") %>%
  filter(selection == "Yes")
head(mydata_filter)


write.table(mydata_filter,"Eastwest_GlobalFst_SWA_SelectionRegions.txt", sep="\t", row.names=FALSE, quote=FALSE)





# Find Loci Within Windows ------------------------------------------------
## read in SWA input file, sorted
mydata_marker_data = read.table("EastvWest/batch_8_SWA_input_eastwest_globalFst_filtered_sorted.txt", header = TRUE, sep = "\t")

## function
find_markers_in_window(marker_data = mydata_marker_data, window_size= 250000, divisions = 150, output = "EastWest_SWA_Num_Name_Loci_Per_Window.txt")




# Filter Loci within Windows by Selection Region --------------------------
## read back in the output from the function above and rename columns
e_find_markers_output <- read.delim("EastWest_GlobalFst_Num_Name_Loci_Per_Window.txt", sep = "\t", header= FALSE, colClasses = c("character", "numeric", "numeric", "character"))
colnames(e_find_markers_output) <- c("chromosome", "position", "num_markers", "loci_names")
View(e_find_markers_output)

## filter data frame with marker names per window to only include selection regions
e_find_markers_filtered <- filter(e_find_markers_output, position %in% mydata_filter$position)
dim(e_find_markers_filtered)
View(e_find_markers_filtered)
dim(mydata_filter)

## write to output file
write.table(e_find_markers_filtered, "EastWest_GlobalFst_SWA_SelectionRegions_Markers.txt", sep = "\t",
            row.names=FALSE, quote=FALSE)




