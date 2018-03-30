############# Plotting Sliding Window Analysis / Individual Loci #########
#
# Adopted from Charlie Waters' script "just_plot_updated_changed_plots"
# To be run AFTER SlidingWindow_MF.R script
# Prepares for & runs plotting function to place all aligned loci over SWA plots
#
# MF edited for Pacific cod 3/26/2018
#
#########################################################################


# Install Packages --------------------------------------------------------

install.packages("ggplot2")
install.packages("dplyr")
library(ggplot2)
library(dplyr)


# Load data ---------------------------------------------------------------
## set working directory
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/SlidingWindow")


## sliding window analysis output
swa_output = read.table("EastvWest/batch_8_SWA_eastwest_output_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150_FILTERED.txt", header = TRUE, sep = "\t")
dim(swa_output)

swa_output_west = read.table("West/batch_8_final_filtered_west_2reg_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150_FILTERED.txt", header = TRUE, sep = "\t")
dim(swa_output_west)


## sliding window analysis input
swa_input = read.table("EastvWest/batch_8_final_filtered_aligned_SWA_input_eastwest.txt", header=TRUE, sep="\t", colClasses = c("numeric", "character", "numeric", "character", "numeric"))
head(swa_input)
colnames(swa_input) <- c("locus", "fst", "chromosome", "position")


swa_input_west = read.table("West/batch_8_SWA_input_west.txt", header=TRUE, sep="\t", colClasses = c("numeric", "character", "numeric", "character", "numeric"))
head(swa_input_west)
colnames(swa_input_west) <- c("locus", "fst", "chromosome", "position")




# Plot --------------------------------------------------------------------
source("Plot_SlidingWindowAnalysis_Functions.R")


plot_outliers(data = swa_output, outlier_data = swa_input, Nb_bootstrap=100000, Nb_divisions = 150, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="plots/EastWest/SWA_EastWest_all_loci_plot")


plot_outliers(data = swa_output_west, outlier_data = swa_input_west, Nb_bootstrap=100000, Nb_divisions = 150, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="plots/West/SWA_West_all_loci_plot")

