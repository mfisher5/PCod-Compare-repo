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
swa_output = read.table("EastvWest/batch_8_SWA_eastwest_output_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150_FILTERED.txt", header = TRUE, sep = "\t")
dim(swa_output)

swa_output_west = read.table("West/batch_8_final_filtered_west_2reg_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150_FILTERED.txt", header = TRUE, sep = "\t")
dim(swa_output_west)


## alignment file in .bed format
align = read.table("../alignment/batch_8_final_filtered_gadMor2LG_filteredMQ_sorted.bed", header=FALSE, sep="\t", colClasses = c("character", "integer", "integer", "character", "integer", "factor"))
dim(align)
colnames(align) <- c("chromosome", "position", "end", "locus", "MappingQuality", "X6")
head(align)



## outlier file. format should be (locus_name \t fst). locus names must match sam file
fstfile = read.table("outlier_data/batch_8_eastwest_Bayescan_p1K_outliers.txt", header = TRUE, colClasses = c("character", "numeric"))
colnames(fstfile) <- c("locus", "fst")
head(fstfile)

fstfile_west = read.table("outlier_data/batch_8_west_OutFLANK_outliers.txt", header = TRUE, colClasses = c("character", "numeric"))
colnames(fstfile_west) <- c("locus", "fst")
head(fstfile_west)





# Join Outlier Location with FST ------------------------------------------
outliers_alldata <- left_join(fstfile,align,by="locus")
outliers <- select(outliers_alldata, c("locus", "fst", "chromosome", "position"))
head(outliers)


outliers_alldata_west <- left_join(fstfile_west,align,by="locus")
outliers_west <- select(outliers_alldata_west, c("locus", "fst", "chromosome", "position"))
head(outliers_west)




# Plot Markers per Window over SWA ----------------------------------------
source("Plot_SlidingWindowAnalysis_Functions.R")
plot_outliers(data = swa_output, outlier_data = outliers, Nb_bootstrap=100000, Nb_divisions = 150, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="SWA_EastWest_outliers_plot")



plot_outliers(data = swa_output_west, outlier_data = outliers_west, Nb_bootstrap=100000, Nb_divisions = 150, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="plots/SWA_West_outliers_plot")
