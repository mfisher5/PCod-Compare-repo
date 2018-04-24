##################### Histograms ##########################
#
# Histograms of Fstats per locus 
#
# MF 4/11/2018
#
###########################################################

library(ggplot2)
library(readxl)
library(plyr)


# Read in Data ------------------------------------------------------------
fstats <- read_excel("D:/Pacific cod/DataAnalysis/PCod-Korea-repo/results/verif/Fstats_perLocus.xlsx")
head(fstats)




# Plot histograms ---------------------------------------------------------

## Fis
ggplot(fstats, aes(x=Fis)) +
  geom_histogram() +
  ylab("Number of Loci") +
  scale_y_continuous(expand=c(0,0), limits = c(0,1600)) +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size=14)) +
  geom_vline(data=fstats, aes(xintercept=mean(fstats$Fis)), linetype="dashed")


## Fst
ggplot(fstats, aes(x=Fst)) +
  geom_histogram() +
  ylab("Number of Loci") +
  scale_y_continuous(expand=c(0,0), limits = c(0,1600)) +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size=14)) +
  geom_vline(data=fstats, aes(xintercept=mean(fstats$Fst)), linetype="dashed")


## Fit
ggplot(fstats, aes(x=Fit)) +
  geom_histogram() +
  ylab("Number of Loci") +
  scale_y_continuous(expand=c(0,0), limits = c(0,1600)) +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size=14)) +
  geom_vline(data=fstats, aes(xintercept=mean(fstats$Fit)), linetype="dashed")








