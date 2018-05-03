# Some background information on AMOVA.
# Analysis of Molecular Variance (AMOVA) is a method of estimating population differentiation 
# directly from molecular data and testing hypotheses about such differentiation. 
# A variety of molecular data - molecular marker data (for example, RFLP or AFLP), 
# direct sequence data, or phylogenetic trees based on such molecular data - 
# may be analyzed using this method (Excoffier, et al. 1992)

# 1- Null hypothesis; samples are taken from a single panmictic population. 
#    Or there is no difference in allele frequencies among the populations.

install.packages("poppr")
library(poppr)
library(pegas)

#The implementation of AMOVA in poppr requires two very basic components: 
# (1) A distance matrix derived from the data and 
# (2) a separate table used to partition the data into different stratifications.
# The distance matrix can be calculated using any distance as long as it is euclidean.

#read in genetic data and hierarchical strata files. 
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/R")

data_all_loci <-read.genepop("../../stacks_b8_wgenome_r05/batch_8_final_filtered_aligned_genepop_eastwest.gen")
levels(pop(data_all_loci))

# Read in a dataframe containing information about the hierarchial levels
my_strata<- read.delim("../AMOVA/EastWest_Strata_for_AMOVA.txt")
head(my_strata)

# assign those levels as strata to the genind object
strata(data_all_loci) <- my_strata

# Change the format of your data from a genind to a genclone object
my_data <- as.genclone(data_all_loci)
my_data

# View some hierarchical levels in the data set
table(strata(my_data, ~population))
table(strata(my_data, ~region, combine = FALSE))
table(strata(my_data, ~site))

# Because the molecular data consist of Euclidean distances derived from vectors of 1s and 0s, 
# the data are unlikely to follow a normal distribution. A null distribution is therefore 
# computed by resampling of the data (Excoffier, et al. 1992). In each permutation, each 
# individual is assigned to a randomly chosen population while holding the sample sizes constant. 
# These permutations are repeated many times, eventually building a null distribution. 
# Hypothesis testing is carried out relative to these resampling distributions.
# amova as implemented in pegas, performs a  hierarchical analysis of molecular variance as described in Excoffier et al. (1992). 
#This implementation accepts any number of hierarchical levels.

#The formula must be of the form d, ~ A/B/... where d is a genclone or genind object, and A, B, etc, 
# are the hierarchical levels from the highest to the lowest one. 
# Any number of levels is accepted, so specifying d ~ A will simply test for population differentiation.
# the poppr.amova function is a wrapper script for amova. It calculates a pairwise distance matrix by default.



###PEGAS AMOVAS#######################################################
#AMOVA considering regions
amova1 <- poppr.amova(my_data, ~region, within = FALSE, method = "pegas", nperm = 1000)
amova1

#The formula must be of the form d, ~ A/B/... where d is a genclone or genind object, 
# and A, B, etc, # are the hierarchical levels from the highest to the lowest one. 
# Any number of levels is accepted, so specifying d ~ A will simply test for population differentiation.
# the poppr.amova function is a wrapper script for amova. It calculates a pairwise distance matrix by default. 


#AMOVA considering regions and populations
amova2<- poppr.amova(my_data, ~region/population, within = FALSE, method = "pegas", 
                            nperm = 1000, quiet = TRUE)

amova2


#### ADE4 AMOVAS##################################################
#AMOVA considering management units

amova3<- poppr.amova(my_data, ~management_unit, within = FALSE, 
                     nperm = 1000, quiet = TRUE)

amova3

set.seed(1999)
amova3.test   <- randtest(amova3, nrepet = 999)


plot(amova3.test)

amova3.test


#AMOVA considering management units and populations


amova4<- poppr.amova(my_data, ~management_unit/population, within = FALSE, 
                     nperm = 1000, quiet = TRUE)

amova4

# To test if populations are significantly different, we perform a randomization test using the function randtest() from the ade4 package. 
#This will randomly permute the sample matrices as described in (Excoffier et al., 1992).

set.seed(1999)
amova4.test   <- randtest(amova4, nrepet = 999)


plot(amova4.test)

amova4.test

#

