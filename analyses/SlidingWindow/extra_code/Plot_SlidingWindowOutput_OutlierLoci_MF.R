############# Plotting Sliding Window Analysis / Outlier Loci #########
#
# Adopted from Charlie Waters' script "just_plot_updated_changed_plots"
# To be run AFTER SlidingWindow_MF.R script
# Prepares for & runs plotting function to place outliers over SWA plots
#
# MF edited for Pacific cod 3/16/2018
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
swa_output = read.table("EastvWest/batch_8_SWA_eastwest_globalFst_filtered_output_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150_FILTERED.txt", header = TRUE, sep = "\t")
dim(swa_output)

swa_output_west = read.table("West/batch_8_final_filtered_west_2reg_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150_FILTERED.txt", header = TRUE, sep = "\t")
dim(swa_output_west)

swa_output_east = read.table("East/batch_8_final_filtered_east_globalFst_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150_FILTERED.txt", header = TRUE, sep = "\t")
dim(swa_output_east)


## SWA input file to get position, fst of outliers
swa_input <- read.table("EastvWest/batch_8_SWA_input_eastwest_globalFst_filtered_sorted.txt", header = TRUE, sep="\t", colClasses = c("character", "numeric", "character", "numeric"))
swa_input_west <- read.table("West/batch_8_SWA_input_west_sorted.txt", header = TRUE, sep="\t", colClasses = c("character", "numeric", "character", "numeric"))
swa_input_east <- read.table("East/batch_8_SWA_input_east_sorted.txt", header = TRUE, sep="\t", colClasses = c("character", "numeric", "character", "numeric"))


## outlier file. format should be (locus_name \t fst). locus names must match sam file
outlierfile = read.table("../outliers/batch_8_final_filtered_aligned_EASTWEST_outliers.txt", header = TRUE, colClasses = c("character", "character"))
outliers <- select(outlierfile, Locus)
colnames(outliers) <- c("Locus")
head(outliers)

outlierfile_west = read.table("../outliers/batch_8_final_filtered_aligned_WEST_outliers.txt", header = TRUE, colClasses = c("character", "numeric"))
outliers_west <- select(outlierfile_west, Locus)
colnames(outliers_west) <- c("Locus")
head(outliers_west)

outlierfile_east = read.table("../outliers/batch_8_final_filtered_aligned_EAST_outliers.txt", header=TRUE, colClasses=c("character", "numeric"))
outliers_east <- select(outlierfile_east, Locus)
colnames(outliers_east) <- c("Locus")
head(outliers_east)




# Join Outlier with Location and FST ------------------------------------------
outdata <- left_join(outliers, swa_input,by="Locus")
head(outdata)

outdata_west <- left_join(outliers_west, swa_input_west,by="Locus")
head(outdata_west)

outdata_east <- left_join(outliers_east, swa_input_east,by="Locus")
head(outdata_west)



# Plot Markers per Window over SWA ----------------------------------------
source("Plot_SlidingWindowAnalysis_Functions.R")
plot_outliers(data = swa_output, outlier_data = outdata, Nb_bootstrap=100000, Nb_divisions = 150, color="black", which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="plots/EastWest/EastWest_SWA_outliers_plot")



plot_outliers(data = swa_output_west, outlier_data = outdata_west, Nb_bootstrap=100000, Nb_divisions = 150, color="deepskyblue4", which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="plots/West/West_SWA_outliers_plot")


plot_outliers(data = swa_output_east, outlier_data = outdata_east, Nb_bootstrap=100000, Nb_divisions = 150, color="mediumorchid2", which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="plots/East/East_SWA_outliers_plot")









# Outliers per Linkage Group ----------------------------------------------







