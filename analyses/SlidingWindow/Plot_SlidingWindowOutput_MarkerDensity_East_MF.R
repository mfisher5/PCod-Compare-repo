############# Plotting Sliding Window Analysis / Marker Density #########
#
# Adopted from Charlie Waters' script "just_plot_updated_changed_plots"
# To be run AFTER SlidingWindow_MF.R script
# Allows for the use of plotting function across ALL chromosomes as input
#
# MF edited for Pacific cod 2/26/2018
#
#########################################################################


# Install Packages --------------------------------------------------------

install.packages("ggplot2")
library(ggplot2)
library(dplyr)

# Load data ---------------------------------------------------------------
## set working directory
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/SlidingWindow")

## input data is the same as for sliding window analysis (the sorted alignment table output as a file in previous script)
swa_input = read.table("batch_8_SWA_input_east_sorted.txt", header = TRUE, sep = "\t")
head(swa_input)



# Run Function to find Marker Density ----------------------------------------------------------
source("SlidingWindow_FindMarkers_Function.R")
find_markers_in_window(marker_data = swa_input, window_size= 250000, divisions = 150, output = "East_Num_Name_Loci_Per_Window.txt")



# Load Marker Density Data from Function Output ---------------------------
## load and subset data file to include only number of markers per window
all_marker_data <- read.table("East_Num_Name_Loci_Per_Window.txt", header=FALSE, sep="\t")
head(all_marker_data)
marker_data <- select(all_marker_data, c("V1", "V2", "V3"))
colnames(marker_data) <- c("chrom", "position", "num_markers")

## create vector of number of markers per window
markers_per_window <- marker_data$num_markers





# Explore Number Markers per Window ------------------------------------------

## find average number of markers per window
mean(markers_per_window) ## 10.97 / 8.11 / 7.17


## plot number markers per window
qplot(markers_per_window, geom="histogram",
      binwidth = 1,
      main = "Markers per Sliding Window, East\nAcross All Linkage Groups",
      xlab = "# Markers", 
      ylab = "# Windows") +
  annotate("text", x = 15, y = 300, label="Mean:\n7.17")




# Plot Markers per Window over SWA ----------------------------------------
## Load in output file from Kot's function to calculate sliding window
swa_output = read.table("batch_8_final_filtered_east_globalFST_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150_FILTERED.txt", header = TRUE, sep = "\t")
dim(swa_output)


## READ IN CHARLIE'S PLOTTING FUNCTION (adjusted)
source("Plot_SlidingWindowAnalysis_Functions.R")
# options: just_plot_swa; just_plot_md; plot_all; plot_outliers; just_plot_overlay



## Run plotting function on all chromosomes
##-- double check dimensions
if(length(swa_output$position) != length(markers_per_window)){
  print("ERROR: position and marker density vectors are not of equal lengths. you cannot run the just-plot function until you subset markers_per_window.")
} else{print("continue to plotting function.")}

##-- plot
just_plot_md(data = swa_output, marker_density = marker_data, Nb_bootstrap=100000, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="plots/East/batch_8_east_SWA_filtered")



## Run plotting function to overlay all loci

plot_outliers(data = swa_output, outlier_data = swa_input, Nb_bootstrap=100000, Nb_divisions = 150, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="plots/East/East_SWA_output_all_loci")









## Run plotting function on specific chromosomes

##-- subset according to # divisions per chromosome
markers_per_window_lg01 <- markers_per_window[1:150]
##-- check dimensions
if(length(sla_output$position) != length(markers_per_window_lg01)){
  print("ERROR: position and marker density vectors are not of equal lengths. you cannot run the just-plot function until you subset markers_per_window.")
} else{print("continue to plotting function.")}
##-- plot
just_plot(data = sla_output, marker_density = output, Nb_bootstrap=100000, which.chromosome.analysis="all", which.chromosome.plot="all",name="batch_8_eastwest_output")
?par









