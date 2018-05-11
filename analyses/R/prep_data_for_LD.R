#################### Prepare data for LD ################
#
# This script will prepare data for use in LD analysis
#   from the `genetics` R package
#
# MF 5/10/2018
#
#########################################################


# Load packages -----------------------------------------------------------
install.packages("adegenet", dependencies = TRUE)
library(adegenet)
library(dplyr)


# Read in data ------------------------------------------------------------
## genepop file


## dataframe of loci and the chromosomes they aligned to
chromosomes <- read.delim("../SlidingWindow/West/batch_8_SWA_input_west_sorted.txt", sep = "\t", header=TRUE)




# Function to write genind object as genepop ------------------------------
## This code was downloaded from https://rdrr.io/github/romunov/zvau/src/R/writeGenPop.R
source(writeGenPop.R)






# Subset genepop by chromosome --------------------------------------------
for(i in unique(chromosomes$LG)){
  chrdata <- filter(chromosomes, chromosome = )
}

