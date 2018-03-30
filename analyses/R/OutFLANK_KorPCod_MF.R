###################### USING OUTFLANK TO IDENTIFY OUTLIERS ##########################

## MF 10/5/2017
## for Korea PCod
## see github: https://github.com/whitlock/OutFLANK/blob/master/R/OutFLANK.R
## PDF with instructions: https://github.com/whitlock/OutFLANK/blob/master/OutFLANK%20readme.pdf


####### Set working directory #######
setwd("D:/Pacific cod/DataAnalysis/PCod-Korea-repo/analyses/Outliers")


####### Install necessary packages ##########

install.packages("devtools")
library(devtools)
source("http://bioconductor.org/biocLite.R")
biocLite("qvalue")
install_github("whitlock/OutFLANK", force = TRUE)
library(OutFLANK)
source("https://raw.githubusercontent.com/whitlock/OutFLANK/master/R/Fst%20Diploids.R")

####### Load input file #######
# File format: input file is a dataframe. Each row is a locus.
#       The column headers are- $LocusName, $FST, $T1, $T2, $FSTNoCorr, $T1NoCorr, $T2NoCorr, $He
# MakeDiploidFSTMat(SNPmat, locusNames, popNames)
#       This will generate the input data frame from an array with rows - individuals, columns - loci. Alleles are coded as (0,1,2,9).

# To convert a genepop file to the input SNPmat, see python script here:https://github.com/mfisher5/PCod-Korea-repo/blob/master/analyses/Outliers/convert_genepop_to_SNPmat.py

loci <- read.table("batch_8_SNPmat_locusnames.txt", header=F)
pops <- read.table("batch_8_SNPmat_popnames.txt", header=F)
data = read.csv("batch_8_final_filtered_SNPmat.txt", header = FALSE, sep = "\t")
View(data)
datamat = as.matrix(data)

FstDataFrame <- MakeDiploidFSTMat(SNPmat = data, locusNames = loci, popNames = pops)





####### Identify Outliers #######


southKOR <- OutFLANK(FstDataFrame, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=9, qthreshold=0.05)

# take a look at the output! 
head(allKOR)
typeof(allKOR) # spoiler: it's a list.
dim(allKOR$results[which(allKOR$results$OutlierFlag==T),]) # how many outliers do you have? 1st #




###### Write output to a text file ########

head(allKOR$results)
outlier_indices <- which(allKOR$results$OutlierFlag == "TRUE")

locus <- c()
he <- c()
fst <-c ()
meanAlleleFreq <- c()
qvals <- c()
pv <- c()
outlier <- c()


for(i in outlier_indices){
  locus <- c(locus, allKOR$results$LocusName[i])
  he <- c(he, allKOR$results$He[i])
  fst <- c(fst, allKOR$results$FST[i])
  meanAlleleFreq <- c(meanAlleleFreq, allKOR$results$meanAlleleFreq[i])
  qvals <- c(qvals, allKOR$results$qvalues[i])
  pv <- c(pv, allKOR$results$pvalues[i])
  outlier <- c(outlier, allKOR$results$OutlierFlag[i])
}

#it's not great, but it works
dataframe <- cbind(locus, he, fst, meanAlleleFreq, qvals, pv, outlier)
write.csv(dataframe, file="allKOR_b8_outflank_outlierstxt")




###### Plotting results #########

# plots the actual (yellow) and theoretical (smoothed blue curve) distribution of Fst
OutFLANKResultsPlotter(allKOR)

# plots Fst against expected Heterozygosity
plot(allKOR$results$FST, allKOR$results$He, xlab = "Per Locus FST", ylab = "Per Locus He")


###########################################################################

loci <- read.table("batch_8_SNPmat_south_locusnames.txt", header=F)
pops <- read.table("batch_8_SNPmat_south_popnames.txt", header=F)
data = read.csv("batch_8_final_filtered_south_SNPmat.txt", header = FALSE, sep = "\t")
View(data)
datamat = as.matrix(data)

FstDataFrame <- MakeDiploidFSTMat(SNPmat = data, locusNames = loci, popNames = PopNames)





####### Identify Outliers #######


southKOR <- OutFLANK(FstDataFrame, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=9, qthreshold=0.05)

# take a look at the output! 
head(southKOR)
typeof(southKOR) # spoiler: it's a list.
dim(southKOR$results[which(southKOR$results$OutlierFlag==T),]) # how many outliers do you have? 1st #




###### Write output to a text file ########

head(southKOR$results)
outlier_indices <- which(southKOR$results$OutlierFlag == "TRUE")

locus <- c()
he <- c()
fst <-c ()
meanAlleleFreq <- c()
qvals <- c()
pv <- c()
outlier <- c()


for(i in outlier_indices){
  locus <- c(locus, southKOR$results$LocusName[i])
  he <- c(he, southKOR$results$He[i])
  fst <- c(fst, southKOR$results$FST[i])
  meanAlleleFreq <- c(meanAlleleFreq, southKOR$results$meanAlleleFreq[i])
  qvals <- c(qvals, southKOR$results$qvalues[i])
  pv <- c(pv, southKOR$results$pvalues[i])
  outlier <- c(outlier, southKOR$results$OutlierFlag[i])
}

#it's not great, but it works
dataframe <- cbind(locus, he, fst, meanAlleleFreq, qvals, pv, outlier)
write.csv(dataframe, file="southKOR_b8_outflank_outliers.txt")




###### Plotting results #########

# plots the actual (yellow) and theoretical (smoothed blue curve) distribution of Fst
OutFLANKResultsPlotter(southKOR)

# plots Fst against expected Heterozygosity
plot(southKOR$results$FST, southKOR$results$He, xlab = "Per Locus FST", ylab = "Per Locus He")



###########################################################################

loci <- read.table("batch_8_SNPmat_EW_locusnames.txt", header=F)
pops <- read.table("batch_8_SNPmat_EW_popnames.txt", header=F)
data = read.csv("batch_8_final_filtered_EW_SNPmat.txt", header = FALSE, sep = "\t")
View(data)
datamat = as.matrix(data)

FstDataFrame <- MakeDiploidFSTMat(SNPmat = data, locusNames = loci, popNames = pops)





####### Identify Outliers #######


ewKOR <- OutFLANK(FstDataFrame, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=9, qthreshold=0.05)

# take a look at the output! 
head(ewKOR)
typeof(ewKOR) # spoiler: it's a list.
dim(ewKOR$results[which(ewKOR$results$OutlierFlag==T),]) # how many outliers do you have? 1st #



###########################################################################

loci <- read.table("batch_8_SNPmat_SW_locusnames.txt", header=F)
pops <- read.table("batch_8_SNPmat_SW_popnames.txt", header=F)
data = read.csv("batch_8_final_filtered_SW_SNPmat.txt", header = FALSE, sep = "\t")
View(data)
datamat = as.matrix(data)

FstDataFrame <- MakeDiploidFSTMat(SNPmat = data, locusNames = loci, popNames = pops)





####### Identify Outliers #######


swKOR <- OutFLANK(FstDataFrame, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=9, qthreshold=0.05)

# take a look at the output! 
head(ewKOR)
typeof(ewKOR) # spoiler: it's a list.
dim(swKOR$results[which(ewKOR$results$OutlierFlag==T),]) # how many outliers do you have? 1st #

