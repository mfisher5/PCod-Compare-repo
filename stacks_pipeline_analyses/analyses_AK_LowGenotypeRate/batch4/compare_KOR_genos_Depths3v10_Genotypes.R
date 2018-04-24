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
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/stacks_pipeline_analyses/analyses_AK_LowGenotypeRate/batch4")

compare_KOR_genos_Depths3v10 <- read_delim("LowStackDepthSim/compare_KOR_genos_Depths3v10_byLocus.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
View(compare_KOR_genos_Depths3v10)

mydata <- compare_KOR_genos_Depths3v10



# Create Mismatch Variable ------------------------------------------------
mydata_plus <- mutate(mydata, mismatched=1-matched)
View(mydata_plus)

mydata_plus_filtered <- mydata_plus[-c(which(mydata_plus$mismatched == 1)),]
max(mydata_plus_filtered$mismatched)

mydata_mismatches <- select(mydata_plus, c(mismatched, diff_homhet, diff_hethom))



# Boxplot -----------------------------------------------------------------

boxplot(mydata_mismatches,  col="darkcyan", xaxt="n", ylab = "Proportion of Shared Loci Genotyped in both Samples", main = "Cause of Differing Genotypes at Locus\n When both Original & Subset Samples Genotyped")
axis(side=1, at=1:3, c("Total Mismatched\n", "\nHeterozygote \nto Homozygote", "\nHomozygote \nto Heterozygote"), tick=FALSE)

boxplot(mydata_mismatches,  col="darkcyan", xaxt="n", ylim = c(0, 0.05), ylab = "Proportion of Shared Loci Genotyped in both Samples", main = "Cause of Differing Genotypes at Locus\n When both Original & Subset Samples Genotyped")
axis(side=1, at=1:3, c("Total Mismatched\n", "\nHeterozygote \nto Homozygote", "\nHomozygote \nto Heterozygote"), tick=FALSE)

mean(mydata_mismatches$mismatched)
# 0.009330867


# Boxplot Missing ---------------------------------------------------------
mydata_missing <- select(mydata, c(n_pairs, both_miss, b1_miss, b2_miss))

boxplot(mydata_missing,  col="darkcyan", xaxt="n", ylab = "Proportion of Shared Loci", main = "Missing Genotypes at Locus")
axis(side=1, at=1:3, c("Missing\nin Both", "Missing\nin Original","Missing \nin Subset"), tick=FALSE)





###################################### RECALLED GENOTYPES ##############################################

# Read in Data ------------------------------------------------------------
compare_KOR_genos_Depths3v10_MB <- read.delim("LowStackDepthSim_revised/KOR_genotypes_depth3v10_MB_perlocus_mismatches_rerun.txt", "\t",header=TRUE)

mydata <- compare_KOR_genos_Depths3v10_MB
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

# Boxplot Mismatches -----------------------------------------------------------------

boxplot(mydata_mismatches,  col="darkcyan", xaxt="n", ylab = "Proportion of Shared Loci Genotyped in both Samples", main = "Cause of Differing Genotypes at Locus\n When both Original & Subset Samples Genotyped")
axis(side=1, at=1:4, c("Total Mismatched\n", "\nHomozygote\nto Homozygote","\nHeterozygote \nto Homozygote", "\nHomozygote \nto Heterozygote"), tick=FALSE)
median(mydata_mismatches$p.mismatched)
mean(mydata_plus$mismatched)




# Boxplot Missing ---------------------------------------------------------
mydata_missing <- select(mydata_plus, c(n.pairs, n.both.miss, n.b1.miss, n.b2.miss))
mydata_missing <- mydata_missing %>%
  mutate(p.both.miss = n.both.miss / n.pairs) %>%
  mutate(p.b1.miss = n.b1.miss / n.pairs) %>%
  mutate(p.b2.miss = n.b2.miss / n.pairs)

mydata_missing <- select(mydata_missing, c(p.both.miss, p.b1.miss, p.b2.miss))

boxplot(mydata_missing,  col="darkcyan", xaxt="n", ylab = "Proportion of Shared Loci", main = "Missing Genotypes at Locus")
axis(side=1, at=1:3, c("Missing\nin Both", "Missing\nin Original","Missing \nin Subset"), tick=FALSE)


