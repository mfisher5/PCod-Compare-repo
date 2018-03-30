
# Load Data ---------------------------------------------------------------

## set working directory
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/SlidingWindow")

## read in files with loci positions and fst
eastwest <- read.delim("EastvWest/batch_8_final_filtered_aligned_SWA_input_eastwest_sorted.txt",header=TRUE)
head(eastwest)
dim(eastwest)

west <- read.delim("West/batch_8_SWA_input_west_sorted.txt",sep="\t", header=TRUE)
head(west)
dim(west)


library(dplyr)

combo <- inner_join(eastwest,west,by="Locus")
head(combo)


plot(combo$fst.x, combo$fst.y, xlab="East v. West Fst", ylab="Within West Fst", main = "Correlation between Within-West and East v. West Fst")
