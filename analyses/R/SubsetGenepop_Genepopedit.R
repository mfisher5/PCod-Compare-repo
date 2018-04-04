############### Subset Genepop with GenepopEdit #################
#
# MF 4/4/2018
# PCod-Compare
#
####################################################################


# Load Packages -----------------------------------------------------------
if(!require("devtools")) install.packages("devtools")
devtools::install_github("rystanley/genepopedit")
install.packages("adegenet")
install.packages("hierfstat")
install.packages("gplots")
install.packages("deldir")
library(deldir)
library(genepopedit)
library(adegenet)
library(hierfstat)
library (gplots)
library(dplyr)


# Import Data -------------------------------------------------------------
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/PCA")

GenePopData <- "batch_8_final_filtered_aligned_genepop_eastwest_GenepopEdit.txt"

## check to see if data file works for genepop edit

metadata <- genepop_detective(GenePopData, variable="All")
metadata$Pops
metadata$PopNum
metadata$Loci
length(metadata$Loci)




# Subset Data -------------------------------------------------------------
## load in outliers file used for Manhattan plots. outlier locus should be first column; other columns don't matter 
outlierfile <- read.table("../outliers/batch_8_final_filtered_aligned_EASTWEST_outliers_snps.txt", sep="\t",header=TRUE,
                          colClasses = c("character", "factor")) # need to change colClasses based on your data columns
head(outlierfile)

## create two lists of outliers
all_outliers <- outlierfile$Locus


## exclude outliers
subset_genepop(genepop= GenePopData, keep = FALSE, subs = all_outliers, path = "batch_8_final_filtered_aligned_genepop_eastwest_neutral.txt")

