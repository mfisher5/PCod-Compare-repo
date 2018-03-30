########################### BAYESCAN OUTPUT SCRIPT ######################
#
# MF 3/9/2018
#
# This script has two options: run the original bayescan plotting funtion, or run an alternative function that plots with ggplot2
#
#
# Bayescan script: will read in Bayescan output files and (1) plot FST v. log10(q), (2) print outlier loci names **will not be same names as your genepop file**, and (3) plot posterior fst distribution
#
# alt. function: must have run the script https://github.com/mfisher5/PCod-Korea-repo/blob/master/analyses/Outliers/bayescan_to_stacks_locus_IDs.py 
#
#
#
###########################################################################


# Install packages --------------------------------------------------------

install.packages("ggplot2")
install.packages("dplyr")
library(ggplot2)
library(dplyr)



# Set Variables -----------------------------------------------------------
#directory with bayescan plotting script (Bayescan script)
bscan_direct <- "D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/R" 

#directory with bayescan output (both)
workingdirect <- "D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/outliers/Bayescan"

#fst file (Bayescan script)
fst_file <- "batch_8_eastwest_p100_output_fst.txt"

#.sel file (Bayescan script)
sel_file <- "batch_8_eastwest_p100_output.sel"

# output file (both)
outfile <- "batch_8_eastwest_p100_outliers.csv"

#set false discovery rate cutoff (both)
fdr <- 0.05

#fst file edited to include stacks locus IDs (alt. function)
fst_edit <- "batch_8_eastwest_Bayescan_p100_fst_stacksIDs.txt"

#plot title (alt. function)
plt.title = "Bayescan Outliers - East v. West (Prior100)"




# Parse / Plot Output with Bayescan Plot Function --------------------------------------------
## my r script is in this directory
setwd(bscan_direct)

## source r script with plot function
source("BAYESCAN_plot_R.r")

## set wd to output file location
setwd(workingdirect)

## produce a list of outlier loci, and a plot of fst v. log10(qvalue)
plot_bayescan(fst_file, FDR = fdr) #FDR sets cutoff

## produce a plot of posterior distribution of any parameter
mydata = read.table(sel_file, colClasses = "numeric")
parameter = "Fst1"
plot(density(mydata[[parameter]]), xlab=parameter, main=paste(parameter,"posterior distribution"))



# ALTERNATIVE: parse & plot with ggplot ---------------------------------------
po <- log10(fdr)

## set wd to output file location
setwd(workingdirect)

## read in fst file
library(readr)
mydata <- read_delim(fst_edit," ", escape_double = FALSE, trim_ws = TRUE)
## edit infile (because of weird spacing in header lines, from bayescan)
head(mydata)

## add column with log10(q val)
mydata <- mutate(mydata, log10q = log10(mydata$qval))
head(mydata)

## add column to identify if loci are outliers or not
mydata <- mutate(mydata, outlier = ifelse(mydata$log10q < po, "True", "False"))
head(mydata)

## output fst file with only outlier loci
outlier_sub <- subset(mydata, outlier == "True")
dim(outlier_sub)
write.csv(outlier_sub, file=outfile, quote=FALSE,  row.names=FALSE)

## plot with ggplot
outliers <- as.factor(mydata$outlier)
base_plot <- ggplot(mydata, aes(x=log10q, y=fst, col=outliers)) +
  geom_point(size=2) +
  scale_x_reverse(lim=c(0.05,-5)) +
  geom_vline(xintercept=po) +
  labs(title=plt.title,y="Fst",x="Log10q") +
  guides(color="none")
base_plot


## plot with ggplot, with labels
labels <- c()
for(i in seq(1, length(mydata$outlier))){
  if(mydata$outlier[i] == "True"){
    labels <- c(labels, mydata$locus[i])
  } else{labels <- c(labels,"")}
}

base_plot +
  geom_text(aes(label=labels), hjust = 1.1, vjust = -.1)

