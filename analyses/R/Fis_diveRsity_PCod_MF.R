# Purpose: getting Fis per population, Fit per locus
#
#From: Mary 10/27/2017
#
#
#
# ---------------------------------------------------------------------------------

# Import these libraries
install.packages("diveRsity")
install.packages("adegenet")
library("diveRsity")
library("adegenet")
library(dplyr)


# Set working directory
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses")


# diveRsity ---------------------------------------------------------------
# Calculate basic descriptive population parameters from a genepop genotype file: Fis per population
mydata_stats <- basicStats(infile = "../stacks_b8_wgenome_r05/batch_8_final_filtered_aligned_genepop_west_byreg.gen")
str(mydata_stats)




# Heterozygosity ---------------------------------------------------------
mydata_het <- mydata_stats$obs_het

pops <- colnames(mydata_het)[2:length(colnames(mydata_het))]

colMeans(mydata_het[,2:length(colnames(mydata_het))])


### If you run option 5 >> suboption 2 in genepop (.DIV output file), the `1-Qintra estimate per sample over all loci with at least one diploid individual typed`` should equal colMeans output






# Adegenet ----------------------------------------------------------------

mydata <-read.genepop("../stacks_b8_verif/batch_8_filteredMAF_filteredIndivids30_filteredLoci_filteredHWE_filteredCR.gen")
mydata_summary <- summary(mydata)
str(mydata_summary)





