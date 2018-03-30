# Purpose: getting Fis per population, Fit per locus
#
#From: Mary 10/27/2017
#
#
#
# ---------------------------------------------------------------------------------

# Import these libraries
install.packages("diveRsity")
library("diveRsity")
install.packages("xlsx")
install.packages("sendplot")
install.packages("plotrix")
install.packages("HWxtest")
library("xlsx")
library("sendplot")
library("plotrix")
library("HWxtest")


# Set working directory
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses")

# Calculate basic descriptive population parameters from a genepop genotype file: Fis per population
basicStats(infile = "../stacks_b4_wgenome/batch_4_MB_filteredMAF_filteredLoci50_filteredIndivids_filteredHWE_noreps.gen", outfile = "stacks_b4_MB_final_filtered_diveRsity_basicStats.txt", fis_ci = TRUE, fis_boots = 1000, fis_alpha = 0.05)


# A minimal function for the calculate of Weir & Cockerham's (2984) Fst and Fit from codominant molecular data: Fit per locus



# ---------------------------------------------------------------------------------

# Set working directory
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses")

# Calculate basic descriptive population parameters from a genepop genotype file: Fis per population
basicStats(infile = "../stacks_b4_wgenome/batch_4_filteredMAF_filteredLoci50_filteredIndivids50_filteredHWE_noreps_genepop.txt", outfile = "b4_final_filtered_diveRsity_basicStats.txt", fis_ci = TRUE, fis_boots = 1000, fis_alpha = 0.05)

