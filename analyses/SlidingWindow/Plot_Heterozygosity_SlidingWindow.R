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

swa_het_output_east <- read.delim("Heterozygosity/batch_8_final_filtered_east_HETEROZYGOSITY_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150.txt")
head(swa_het_output_east)

swa_het_output_west <- read.delim("Heterozygosity/batch_8_final_filtered_west_HETEROZYGOSITY_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150.txt")
head(swa_het_output_west)




# Plot just heterozygosity --------------------------------------------------------------------
source("Plot_SlidingWindowAnalysis_Functions.R")

plot_het_overlay(data1 = swa_het_output_east, data2 = swa_het_output_west, legend.text.data1 = "East", legend.text.data2 = "West", Nb_bootstrap=100000, Nb_divisions = 150, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="Heterozygosity/SWA_Het_East_West_overlay")




# Plot Het under East v. West FST -----------------------------------------
source("Plot_SlidingWindowAnalysis_FstHet_Function.R")

swa_output = read.table("EastvWest/batch_8_SWA_eastwest_globalFst_filtered_output_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150_FILTERED.txt", header = TRUE, sep = "\t")

plot_fst_het(data = swa_output,
             het1 = swa_het_output_east,
             het2 = swa_het_output_west,
             legend.text.het1 = "East", legend.text.het2 = "West",
             Nb_bootstrap_fst=100000, Nb_divisions_fst = 150, 
             Nb_bootstrap_het=100000, Nb_divisions_het = 150,
             which.chromosome.analysis="all", which.chromosome.plot="all",
             export = TRUE, name="plots/EastWest/EastvWest_Het_filtered")




