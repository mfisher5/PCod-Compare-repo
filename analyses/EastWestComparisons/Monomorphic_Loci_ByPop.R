




# Load Packages -----------------------------------------------------------
install.packages("adegenet")
library(adegenet)
mydata <- read.genepop(file = "batch_8_final_filtered_aligned_genepop.gen")
basicstat <- summary(mydata)
str(basicstat)



if (!require("devtools")) install.packages("devtools") # to install
#install the package from *Github*
devtools::install_github("rystanley/genepopedit") 
library(genepopedit) # load the library

library(ggplot2)
library(dplyr)



# Set Working Directory ---------------------------------------------------
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/stacks_b8_wgenome_r05")



# Calculate Major Allele Frequency by Pop ----------------------------------------

allele_freqs_bypop <- genepop_allelefreq("batch_8_final_filtered_aligned_genepop_forGenepopEdit.txt")
genepop_detective("batch_8_final_filtered_aligned_genepop_forGenepopEdit.txt", variable="Pops")

ak_pops <- c("KOD03","AD06","WC05","HS04","SS1213","JDF12","PWS12","UP03")
kor_pops <- c("PO2015","GE2015","NA021015","YS121315","JUK07","JB121807","JB021108","BOR07","GEO020414")

maf_ak <- filter(allele_freqs_bypop, Population %in% ak_pops)
hist(maf_ak$MAF)

maf_kor <- filter(allele_freqs_bypop, Population %in% kor_pops)
hist(maf_kor$MAF)

ggplot(maf_ak,aes(x=MAF, fill = Population)) + 
  geom_histogram(data = maf_ak, bins= 20) +
  facet_wrap(~Population) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(maf_kor,aes(x=MAF, fill = Population)) + 
  geom_histogram(data = maf_kor, bins= 20) +
  facet_wrap(~Population) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## filter for maf = 1.0
kor_monomorphic <- maf_kor %>%
  filter(MAF == 1.0) %>%
  count(Population)

ak_monomorphic <- maf_ak %>%
  filter(MAF == 1.0) %>%
  count(Population)
kor_monomorphic
ak_monomorphic






# Major Allele Frequency Across Pops --------------------------------------
maf_ak_grouped <- maf_ak %>%
  group_by(Loci) %>%
  summarise(avg = mean(MAF))

maf_list_ak <- maf_ak_grouped$avg
length(which(maf_list_ak == 1.0))



maf_kor_grouped <- maf_kor %>%
  group_by(Loci) %>%
  summarise(avg = mean(MAF))

maf_list_kor <- maf_kor_grouped$avg
length(which(maf_list_kor == 1.0))




