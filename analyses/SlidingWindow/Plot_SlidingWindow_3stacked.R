############## Plot Full Data, East and West Together ##############
#
# To plot three sliding window analyses for LG19
# 
# MF 4/6/2018 for PCod Compare Project
#
####################################################################


# Install Packages --------------------------------------------------------

install.packages("ggplot2")
library(ggplot2)
library(dplyr)
source("Plot_SlidingWindowAnalysis_Functions.R")

# Load data ---------------------------------------------------------------
## set working directory
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/SlidingWindow")

east <- read.table("East/batch_8_final_filtered_east_globalFST_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150_FILTERED.txt", header=TRUE, sep="\t")
head(east)

west <- read.table("West/batch_8_final_filtered_west_2reg_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150_FILTERED.txt", header=TRUE, sep="\t")
head(west)

all <- read.table("EastvWest/batch_8_SWA_eastwest_globalFst_filtered_output_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150_FILTERED.txt", header=TRUE, sep="\t")




# Plot just Three Lines ---------------------------------------------------

plot_lines_overlay(data1 = all, data2 = east, data3 = west, Nb_divisions = 150, legend.text=c("All Data", "East Population", "West Population"), which.chromosome.analysis="all", which.chromosome.plot="all",export = FALSE, name="SWA_lines_overlay_test")



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


all_selection <- all %>%
  mutate(positive = ifelse(`Fst.Fct` > upper_95, "1", "0")) %>%
  mutate(negative = ifelse(`Fst.Fct` < lower_95, "1", "0")) %>%
  mutate(selection = ifelse(positive == 1, "Yes", "No"))
head(all_selection)




# Plot All Lines with Divergence regions ----------------------------------

plot_lines_overlay_divergence(data1 = all_selection, data2 = east_selection, data3 = west_selection, Nb_divisions = 150, legend.text=c("All Data", "East Population", "West Population"), which.chromosome.analysis="all", which.chromosome.plot="all",export = FALSE, name="plots/SWA_lines_overlay_divergence_test")




