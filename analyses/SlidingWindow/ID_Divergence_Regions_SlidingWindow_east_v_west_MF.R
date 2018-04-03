################# Identify Regions under Selection from Sliding Window Analysis #########
#
# This script will filter the output file from sliding window analysis...
#    to include only regions where the FST average goes above or below the confidence interval
# You can then compare overlapping regions between two analyses on plots per linkage group
# You can also identify which loci are within the selected regions
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

east <- read_delim("East/batch_8_final_filtered_east_globalFST_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150_FILTERED.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
head(east)

west <- read_delim("West/batch_8_final_filtered_west_2reg_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150_FILTERED.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
head(west)



# Add columns for selection -------------------------------------------------------------
east_selection <- east %>%
  mutate(positive = ifelse(`Fst.Fct` > upper_95, "1", "0")) %>%
  mutate(negative = ifelse(`Fst.Fct` < lower_95, "1", "0")) %>%
  mutate(selection = ifelse(positive == 1, "Yes", "No"))
head(east_selection)

west_selection <- west %>%
  mutate(positive = ifelse(`Fst.Fct` > upper_95, "1", "0")) %>%
  mutate(negative = ifelse(`Fst.Fct` < lower_95, "1", "0")) %>%
  mutate(selection = ifelse(positive == 1, "Yes", "No"))
head(west_selection)




# Visualize Fst v. Bootstrap CI ----------------------------------------------
## positive selection: Fst above upper 95% confidence interval
ggplot(data=east_selection, aes(x = `Fst.Fct`, y = upper_95)) +
  geom_point(aes(color=positive)) +
  labs(title="East Sliding Window\nPositive Selection")

## negative selection: Fst below lower 95% confidence interval
ggplot(data=east_selection, aes(x = `Fst.Fct`, y = lower_95)) +
  geom_point(aes(color=negative)) +
  labs(title="East Sliding Window\nNegative Selection")


## positive selection: Fst above upper 95% confidence interval
ggplot(data=west_selection, aes(x = `Fst.Fct`, y = upper_95)) +
  geom_point(aes(color=positive)) +
  ylim(-0.1,0.5) +
  xlim(-0.1,0.5) +
  labs(title="West Sliding Window\nPositive Selection")

## negative selection: Fst below lower 95% confidence interval
ggplot(data=west_selection, aes(x = `Fst.Fct`, y = lower_95)) +
  geom_point(aes(color=negative)) +
  ylim(-0.1,0.1) +
  xlim(-0.1,0.5) +
  labs(title="West Sliding Window\nNegative Selection")








# Visualize P values ------------------------------------------------------
ggplot(data=east_selection, aes(x = `Fst.Fct`, y = pvalue)) +
  geom_point(aes(color=selection)) +
  ylim(-0.1,1) +
  xlim(-0.1,1) +
  labs(title="East Sliding Window\nP values")

ggplot(data=west_selection, aes(x = `Fst.Fct`, y = pvalue)) +
  geom_point(aes(color=selection)) +
  ylim(-0.1,1) +
  xlim(-0.1,0.5) +
  labs(title="West Sliding Window\nP values")






# Match Data sets on Plot ---------------------------------------------------------
## plot all chromosomes using overlay function
colnames(east_selection) <- c("chromosome","position","Fst/Fct","Mean_boostrap","lower_95","upper_95","pvalue","positive","negative","selection")
colnames(west_selection) <- c("chromosome","position","Fst/Fct","Mean_boostrap","lower_95","upper_95","pvalue","positive","negative","selection")

just_plot_overlay_diverge(data1 = east_selection, data2 = west_selection, Nb_divisions = 150, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="plots/SWA_East_West_Divergence_Overlay")







# Filter Data and Write to File --------------------------------------------
east_filter <- east_selection %>%
  filter(selection == "Yes")
dim(east_filter)

west_filter <- west_selection %>%
  filter(selection == "Yes")
dim(west_filter)



# Find Loci Within Windows ------------------------------------------------
## read in SWA input file, sorted
east_marker_data = read.table("batch_8_SWA_input_east_sorted.txt", header = TRUE, sep = "\t")
west_marker_data = read.table("batch_8_SWA_input_west_sorted.txt", header = TRUE, sep = "\t")

## function
find_markers_in_window(marker_data = east_marker_data, window_size= 250000, divisions = 150, output = "east_SWA_Num_Name_Loci_Per_Window.txt")

find_markers_in_window(marker_data = west_marker_data, window_size= 250000, divisions = 150, output = "West_SWA_Num_Name_Loci_Per_Window.txt")




# Filter Loci within Windows by Selection Region --------------------------
## east data
e_find_markers_output <- read.delim("East_Num_Name_Loci_Per_Window.txt", sep = "\t", header= FALSE, colClasses = c("character", "numeric", "numeric", "character"))
colnames(e_find_markers_output) <- c("chromosome", "position", "num_markers", "loci_names")
View(e_find_markers_output)

e_find_markers_filtered <- filter(e_find_markers_output, position %in% east_filter$position)
dim(e_find_markers_filtered)
View(e_find_markers_filtered)
dim(east_filter)

write.table(e_find_markers_filtered, "East_SWA_SelectionRegions_Markers.txt", sep = "\t", 
            row.names=FALSE, quote=FALSE)

## west data
w_find_markers_output <- read.delim("West_SWA_Num_Name_Loci_Per_Window.txt", sep = "\t", header= FALSE, colClasses = c("character", "numeric", "numeric", "character"))
colnames(w_find_markers_output) <- c("chromosome", "position", "num_markers", "loci_names")
View(w_find_markers_output)

w_find_markers_filtered <- filter(w_find_markers_output, position %in% west_filter$position)
dim(w_find_markers_filtered)
dim(west_filter)

write.table(w_find_markers_filtered, "West_SWA_SelectionRegions_Markers.txt", sep = "\t", 
            row.names=FALSE, quote=FALSE)




