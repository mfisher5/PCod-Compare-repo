################# PCadapt ##################
#
# MF 3/30/2018
#
# Reference: https://github.com/bcm-uga/pcadapt/blob/master/vignettes/pcadapt.Rmd
#
#############################################


# Load packages -----------------------------------------------------------

install.packages("devtools")
install.packages("curl")
library(curl)
devtools::install_github("bcm-uga/pcadapt")
source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")
library(pcadapt)
library(qvalue)




# Load in Data & Set File Names -------------------------------------------------------------
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/outliers/PCAdapt")
## VCF file
mydata <- read.pcadapt("batch_8_eastwest.vcf", type = "vcf")

## file with list of loci IDs, same order as in VCF file
myloci <- read.table("batch_8_eastwest_vcf_loci.txt", header=FALSE)

## output file with outlier information
outfile <- "batch_8_eastwest_PCAdapt_q005_k2_outliers.txt"

## alpha cutoff for q value
alpha <- 0.05
## for scoreplot: create list of population names in same order as VCF
#poplist.names <- c(rep("POH", 31), rep("GE", 50), rep("NAM", 11), rep("JinBay", 81), rep("YSBlock", 25), rep("BOR", 21))
poplist.names <- c(rep("Adak", 41), rep("WaCo", 41), rep("HStrait", 46), rep("PWSound", 47), rep("Unimak", 42))


# Choose Principal Components ---------------------------------------------
x <- pcadapt(mydata, K = 20, method = "mahalanobis", min.maf = 0.05)

## with scree plot
plot(x, option="screeplot") #all principal components
plot(x, option = "screeplot", K = 10) #first 10 PCs

## with score plot
plot(x, option="scores", pop=poplist.names) #plot PCA




# Compute Test Statistic --------------------------------------------------
mydata.pcadapt <- pcadapt(mydata, K = 3)

plot(mydata.pcadapt, option="manhattan")
hist(mydata.pcadapt$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(mydata.pcadapt, option = "stat.distribution")





# Use q-values to ID outliers ---------------------------------------------

# make q values from p values
qval <- qvalue(mydata.pcadapt$pvalues)$qvalues

# apply alpha cutoff & return indices of marker numbers
outliers <- which(qval < alpha)
length(outliers)
outliers # THIS IS THE INDEX OF THE LOCUS, NOT THE LOCUS NAME


# Use Benjamini-Hochberg Procedure to ID outliers -------------------------

padj <- p.adjust(mydata.pcadapt$pvalues,method="BH")
outliers <- which(padj < alpha)
length(outliers)
outliers

max(mydata.pcadapt$pvalues, na.rm = TRUE)
mydata.pcadapt$pvalues[335]


# Match Outlier Index to Locus Name ---------------------------------------
myloci$V1 <- as.character(myloci$V1) #important if loci IDs are `locus_snp`
head(myloci)

## if using q values procedure
outlier_IDs <- c()
pvals <- c()
qvals <- c()
dist <- c()
for(i in outliers){
  outlier_IDs <- c(outlier_IDs, myloci$V1[i])
  pvals <- c(pvals, mydata.pcadapt$pvalues[i])
  qvals <- c(qvals, qval[i])
  dist <- c(dist, mydata.pcadapt$stat[i])
}
outlier_data <- cbind(outliers, outlier_IDs, pvals, qvals, dist)
colnames(outlier_data) <- c("VCF_Index", "Locus", "Pvalue", "Qvalue", "sq.M.dist")
View(outlier_data)


write.table(outlier_data, outfile, quote=FALSE,
            sep="\t", row.names = FALSE)




