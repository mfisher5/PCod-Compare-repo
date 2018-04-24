################ plot per pop Fis for Subset Simulation #############
#
# MF 4/23/2018
#
######################################################################

library(ggplot2)
library(dplyr)
library(reshape2)

# Read in Data ------------------------------------------------------------
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/stacks_pipeline_analyses/analyses_AK_LowGenotypeRate/batch4/LowStackDepthSim_revised")

OS <- read.delim("Fis_perpop_batch1_Stacks.txt")
OR <- read.delim("Fis_perpop_batch1_MB.txt")
SS <- read.delim("Fis_perpop_batch2_Stacks.txt")
SR <- read.delim("Fis_perpop_batch2_MB.txt")



# Combine data ---------------------------------------------------------------
OS <- mutate(OS, Data = "Stacks, Original")
OR <- mutate(OR, Data = "Recalled, Original")
SS <- mutate(SS, Data = "Stacks, Subset")
SR <- mutate(SR, Data = "Recalled, Subset")
all_df <- rbind(OS, OR, SS, SR)
View(all_df)



# Plot --------------------------------------------------------------------
ggplot(all_df, aes(x=Population, y=Fis)) +
  geom_col() +
  xlab("Sampling Site") +
  facet_wrap(~Data) +
  theme(axis.title = element_text(size=12, face = "bold"), axis.text.y = element_text(size=12), 
        axis.text.x = element_text(size = 12, angle=90, hjust = 1),
        strip.text= element_text(size = 12, face = "bold"))
  

