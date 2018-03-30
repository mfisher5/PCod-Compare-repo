## This script plots the differences in genotypes between Korean samples (full file size, m = 10) and subset Korean samples (half file size, m = 3) ##

library(readr)
compare_KOR_genos_Depths3v10 <- read_delim("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/stacks_pipeline_analyses/analyses_AK_LowGenotypeRate/compare_KOR_genos_Depths3v10.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
View(compare_KOR_genos_Depths3v10)

data <- compare_KOR_genos_Depths3v10

p_genotypedb1 = data$genotyped.b1 / 11496

p_genotypedb2 = data$genotyped.b2 / 11496

p_diff = data$diff.genos / min(data$genotyped.b1, data$genotyped.b2)

#combine datasets into dataframe
p_combined = data.frame(p_genotypedb1, p_genotypedb2)

#plot
boxplot(p_combined, col="darkcyan", names = c("Genotyped in Orig.", "Genotyped in Subset"), ylab = "Proportion of Shared Loci", main = "Loci that were Genotyped in Original v. Subset Korean data")

#######

p_diff = data$diff.genos / min(data$genotyped.b1, data$genotyped.b2)

#create list of proportion loci different b/c of missing data in original samples
diff_b1miss = c()
for(i in seq(1, length(data$diff.miss_het))){
  new_diff = (data$diff.miss_het[i] + data$diff.miss_hom[i]) / data$diff.genos[i]
  diff_b1miss = c(diff_b1miss, new_diff)
} 
#create list of proportion loci different b/c of missing data in subset samples
diff_b2miss = c()
for(i in seq(1, length(data$diff.het_miss))){
  new_diff = (data$diff.het_miss[i] + data$diff.hom_miss[i]) / data$diff.genos[i]
  diff_b2miss = c(diff_b2miss, new_diff)
} 

p_diff_missing_combined <- data.frame(diff_b1miss, diff_b2miss)
#plot
boxplot(p_diff_missing_combined, col="darkcyan", names = c("Mismatch: Original Genotype Missing", "Mismatch: Subset Genotype Missing"), ylab = "Proportion of Loci Genotyped that Did Not Match", main = "Differing Genotypes at Locus\n From Missing Data")

###### 

#create list of proportions of genotypes different b/c orig. heterozygote, subset homozygote
p_hethom = c()
for(i in seq(1, length(data$diff.het_hom))){
  new_diff = data$diff.het_hom[i] / min(data$genotyped.b1[i], data$genotyped.b2[i])
  p_hethom[i] = new_diff
}
#create list of proportions of genotypes different b/c orig. homozygote, subset heterozygote
p_homhet = c()
for(i in seq(1, length(data$diff.hom_het))){
  new_diff = data$diff.hom_het[i] / min(data$genotyped.b1[i], data$genotyped.b2[i])
  p_homhet[i] = new_diff
}

#create list of proportions of genotypes different where both samples were genotyped
p_diff_true <- c()
for(i in seq(1, length(data$diff.hom_het))){
  new_diff = (data$diff.hom_het[i] + data$diff.hom_het[i]) / min(data$genotyped.b1[i], data$genotyped.b2[i])
  p_diff_true[i] = new_diff
}


p_diff_combined <- data.frame(p_diff_true, p_hethom, p_homhet)

#plot
boxplot(p_diff_combined, col="darkcyan", names = c("All Mismatched Genotypes", "Orig. Het -> Subset Hom", "Orig. Hom -> Subset Het"), ylab = "Proportion of Loci Genotyped in both Batches", main = "Cause of Differing Genotypes at Locus\n When both Batches Genotyped")

####################################################################################
library(readr)
compare_KOR_genos_Depths3v10_MBcorrect <- read_delim("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/stacks_pipeline_analyses/analyses_AK_LowGenotypeRate/compare_KOR_genos_Depths3v10_MBcorrect.txt",
"\t", escape_double = FALSE, trim_ws = TRUE)
View(compare_KOR_genos_Depths3v10_MBcorrect)


data <- compare_KOR_genos_Depths3v10_MBcorrect


p_genotypedb1 = data$genotyped.b1 / 6937

p_genotypedb2 = data$genotyped.b2 / 6937

p_diff = data$diff.genos / min(data$genotyped.b1, data$genotyped.b2)

#combine datasets into dataframe
p_combined = data.frame(p_genotypedb1, p_genotypedb2)

#plot
boxplot(p_combined, col="darkcyan", names = c("Genotyped in Orig.", "Genotyped in Subset"), ylab = "Proportion of Shared Loci", main = "Loci that were Genotyped in Original v. Subset Korean data\nMB Corrections")

######


#create list of proportion loci different b/c of missing data in original samples
diff_b1miss = c()
for(i in seq(1, length(data$diff.miss_het))){
  new_diff = (data$diff.miss_het[i] + data$diff.miss_hom[i]) / data$diff.genos[i]
  diff_b1miss = c(diff_b1miss, new_diff)
} 
#create list of proportion loci different b/c of missing data in subset samples
diff_b2miss = c()
for(i in seq(1, length(data$diff.het_miss))){
  new_diff = (data$diff.het_miss[i] + data$diff.hom_miss[i]) / data$diff.genos[i]
  diff_b2miss = c(diff_b2miss, new_diff)
} 

p_diff_missing_combined <- data.frame(diff_b1miss, diff_b2miss)
#plot
boxplot(p_diff_missing_combined, col="darkcyan", names = c("Mismatch: Original Genotype Missing", "Mismatch: Subset Genotype Missing"), ylab = "Proportion of Loci Genotyped that Did Not Match", main = "Differing Genotypes at Locus From Missing Data\nMB Correction")

###### 

#create list of proportions of genotypes different b/c orig. heterozygote, subset homozygote
p_hethom = c()
for(i in seq(1, length(data$diff.het_hom))){
  new_diff = data$diff.het_hom[i] / min(data$genotyped.b1[i], data$genotyped.b2[i])
  p_hethom[i] = new_diff
}
#create list of proportions of genotypes different b/c orig. homozygote, subset heterozygote
p_homhet = c()
for(i in seq(1, length(data$diff.hom_het))){
  new_diff = data$diff.hom_het[i] / min(data$genotyped.b1[i], data$genotyped.b2[i])
  p_homhet[i] = new_diff
}

#create list of proportions of genotypes different where both samples were genotyped
p_diff_true <- c()
for(i in seq(1, length(data$diff.hom_het))){
  new_diff = (data$diff.hom_het[i] + data$diff.hom_het[i] + data$diff.het_het + data$diff.hom_hom) / min(data$genotyped.b1[i], data$genotyped.b2[i])
  p_diff_true[i] = new_diff
}


p_diff_combined <- data.frame(p_diff_true, p_hethom, p_homhet)

#plot
boxplot(p_diff_combined, col="darkcyan", names = c("All Mismatched Genotypes", "Orig. Het -> Subset Hom", "Orig. Hom -> Subset Het"), ylab = "Proportion of Loci Genotyped in both Batches", main = "Cause of Differing Genotypes at Locus When both Batches Genotyped\n (MB Correction)")

####################################################################################