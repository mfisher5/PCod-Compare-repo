########## Use Hierfstat to calculate Ho per Locus #############
#
# MF 3/29/2018
#
###############################################################


install.packages("hierfstat")
install.packages("adegenet")
library(hierfstat)
library(adegenet)
library(stringr)
library(dplyr)

setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo")
mydata <- read.genepop(file = "stacks_b8_wgenome_r05/batch_8_final_filtered_aligned_genepop.gen")

basic_stats_all <- basic.stats(mydata)


# Prep SWA input for west -------------------------------------------------
kordata <- read.genepop(file = "stacks_b8_wgenome_r05/batch_8_final_filtered_aligned_genepop_west.gen")
basic_stats_kor <- basic.stats(kordata)

ho_data <- basic_stats_kor$perloc %>%
  select(Ho)
Locus <- rownames(ho_data)
ho_data <- cbind(Locus, ho_data)
ho_data$Locus <- ho_data$Locus %>% str_replace("_.*","") %>% str_replace("X","")
ho_data$Locus <- as.integer(ho_data$Locus)
View(ho_data)


write.table(ho_data, "analyses/SlidingWindow/West_Ho_perLocus.txt", row.names=FALSE, quote = FALSE)



# Prep SWA input for east -------------------------------------------------
akdata <- read.genepop(file = "stacks_b8_wgenome_r05/batch_8_final_filtered_aligned_genepop_east_coastal.gen")
basic_stats_ak <- basic.stats(akdata)

ho_data <- basic_stats_ak$perloc %>%
  select(Ho)
Locus <- rownames(ho_data)
ho_data <- cbind(Locus, ho_data)
ho_data$Locus <- ho_data$Locus %>% str_replace("_.*","") %>% str_replace("X","")
ho_data$Locus <- as.integer(ho_data$Locus)
View(ho_data)


write.table(ho_data, "analyses/SlidingWindow/East_Ho_perLocus.txt", row.names=FALSE, quote = FALSE)

