###################### USING OUTFLANK TO IDENTIFY OUTLIERS ##########################

## MF 10/5/2017
## for Korea PCod
## see github: https://github.com/whitlock/OutFLANK/blob/master/R/OutFLANK.R
## PDF with instructions: https://github.com/whitlock/OutFLANK/blob/master/OutFLANK%20readme.pdf


####### Set working directory #######
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/Outliers/OutFLANK")


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

loci <- read.table("batch_8_east_SNPmat_locusnames.txt", header=F)
pops <- read.table("batch_8_east_SNPmat_popnames.txt", header=F)
data = read.csv("batch_8_east_SNPmat.txt", header = FALSE, sep = "\t")
datamat = as.matrix(data)

FstDataFrame <- MakeDiploidFSTMat(SNPmat = data, locusNames = loci, popNames = pops)





####### Identify Outliers #######


mydata.outflank <- OutFLANK(FstDataFrame, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples=6, qthreshold=0.05)

# take a look at the output! 
typeof(mydata.outflank) # spoiler: it's a list.
dim(mydata.outflank$results[which(mydata.outflank$results$OutlierFlag==T),]) # how many outliers do you have? 1st #




###### Write output to a text file ########

head(mydata.outflank$results)
outlier_indices <- which(mydata.outflank$results$OutlierFlag == "TRUE")

locus <- c()
he <- c()
fst <-c ()
meanAlleleFreq <- c()
qvals <- c()
pv <- c()
outlier <- c()

# 3/25 - fixed bug where stacks locus IDs with "_" were being written out incorrectly
for(i in outlier_indices){
  locus <- c(locus, as.character(mydata.outflank$results$LocusName[i]))
  print(as.character(mydata.outflank$results$LocusName[i]))
  he <- c(he, mydata.outflank$results$He[i])
  fst <- c(fst, mydata.outflank$results$FST[i])
  meanAlleleFreq <- c(meanAlleleFreq, mydata.outflank$results$meanAlleleFreq[i])
  qvals <- c(qvals, mydata.outflank$results$qvalues[i])
  pv <- c(pv, mydata.outflank$results$pvalues[i])
  outlier <- c(outlier, mydata.outflank$results$OutlierFlag[i])
}

#write out to csv file
dataframe <- cbind(locus, he, fst, meanAlleleFreq, qvals, pv, outlier)
write.csv(dataframe, file="batch_8_east_bypop_outflank_outliers.csv", row.names = FALSE, quote=FALSE)




###### Plotting results #########

# plots the actual (yellow) and theoretical (smoothed blue curve) distribution of Fst
OutFLANKResultsPlotter(mydata.outflank)

# plots Fst against expected Heterozygosity
plot(mydata.outflank$results$FST, mydata.outflank$results$He, xlab = "Per Locus FST", ylab = "Per Locus He")


