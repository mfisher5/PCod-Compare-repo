################## Testing Genotype Error with Lower Stack Depth ################
#
## This script plots the differences in genotypes between Korean samples (full file size, m = 10) and subset Korean samples (half file size, m = 3) 
## Note that the python Jupyter Notebook 'Testing Lower Stack Depth by Locus Mismatches' must be run before this R code to produce the input files
#
# MF, Last updated 4/23/2018
#
#################################################################################

library(readr)
library(ggplot2)
library(dplyr)


###################################### STACKS CALLING GENOTYPES ############################

# Read in Data ------------------------------------------------------------
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/stacks_pipeline_analyses/analyses_AK_LowGenotypeRate/batch4/LowStackDepthSim_revised")

mydata <- read_delim("KOR_genotypes_depth3v10_byIndivid_mismatches.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
View(mydata)

# Create Mismatch Variable ------------------------------------------------
mydata_plus <- mutate(mydata, p.matched = n.matched / n.both.genod)
View(mydata_plus)

mydata_plus <- mydata_plus %>%
  mutate(p.mismatched = 1-p.matched) %>%
  mutate(p.mismatch.homhom = n.mismatch.hom2hom / n.both.genod) %>%
  mutate(p.mismatch.homhet = n.mismatch.hom2het / n.both.genod) %>%
  mutate(p.mismatch.hethom = n.mismatch.het2hom / n.both.genod)

mydata_mismatches <- select(mydata_plus, c(p.mismatched, p.mismatch.homhom, p.mismatch.hethom, p.mismatch.homhet))
mydata_mismatches <- filter(mydata_mismatches, !is.na(p.mismatched))

mean(mydata_mismatches$p.mismatched)



# Boxplot -----------------------------------------------------------------

boxplot(mydata_mismatches,  col="darkcyan", xaxt="n", ylab = "Proportion of Shared Loci Genotyped in both Samples", main = "Cause of Differing Genotypes in Individual\n When both Original & Subset Samples Genotyped")
axis(side=1, at=1:4, c("Total Mismatched\n", "Homozygote\nto Homozygote","\nHeterozygote \nto Homozygote", "\nHomozygote \nto Heterozygote"), tick=FALSE)


mean(mydata_mismatches$p.mismatched)
# 0.003265968






###################################### RECALLED GENOTYPES ##############################################

# Read in Data ------------------------------------------------------------
mydata <- read.delim("KOR_genotypes_depth3v10_MB_byIndivid_mismatches.txt", "\t",header=TRUE)

View(mydata)


# Create Mismatch Variable ------------------------------------------------
mydata_plus <- mutate(mydata, p.matched = n.matched / n.both.genod)
View(mydata_plus)

mydata_plus <- mydata_plus %>%
  mutate(p.mismatched = 1-p.matched) %>%
  mutate(p.mismatch.homhom = n.mismatch.hom2hom / n.both.genod) %>%
  mutate(p.mismatch.homhet = n.mismatch.hom2het / n.both.genod) %>%
  mutate(p.mismatch.hethom = n.mismatch.het2hom / n.both.genod)

mydata_mismatches <- select(mydata_plus, c(p.mismatched, p.mismatch.homhom, p.mismatch.hethom, p.mismatch.homhet))
mydata_mismatches <- filter(mydata_mismatches, !is.na(p.mismatched))

mean(mydata_mismatches$p.mismatched)
## 0.07794633

# Boxplot Mismatches -----------------------------------------------------------------

boxplot(mydata_mismatches,  col="darkcyan", xaxt="n", ylab = "Proportion of Shared Loci Genotyped in both Samples", main = "Cause of Differing Genotypes in Individual\n When both Original & Subset Samples Genotyped")
axis(side=1, at=1:4, c("Total Mismatched\n", "\nHomozygote\nto Homozygote","\nHeterozygote \nto Homozygote", "\nHomozygote \nto Heterozygote"), tick=FALSE)



