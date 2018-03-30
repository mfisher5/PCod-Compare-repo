############## Plot Heterozygosity for East v. West ################
#
# MF 3/30/2018
#
####################################################################




# Install Packages --------------------------------------------------------

install.packages("TTR")
library(TTR) 
library(zoo) # don't need to install; will install with TTR

install.packages("dplyr") #always install this last!
library(dplyr)


# Load Data ---------------------------------------------------------------

## set working directory
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/SlidingWindow")

swa_het_output_east <- read.delim("batch_8_final_filtered_east_HETEROZYGOSITY_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150.txt")
head(swa_het_output_east)

swa_het_output_west <- read.delim("batch_8_final_filtered_west_HETEROZYGOSITY_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150.txt")
head(swa_het_output_west)




# Plot --------------------------------------------------------------------
source("Plot_SlidingWindowAnalysis_Functions.R")

plot_het_overlay(data1 = swa_het_output_east, data2 = swa_het_output_west, legend.text.data1 = "East", legend.text.data2 = "West", Nb_bootstrap=100000, Nb_divisions = 150, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="SWA_Het_East_West_overlay")


