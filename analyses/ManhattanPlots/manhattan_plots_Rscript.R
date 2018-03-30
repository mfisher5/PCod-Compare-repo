################# Manhattan Plot Code #################
#
# Modified from Charlie Waters' script
#
# MF 2/8/2018
#
#######################################################

############################  BATCH 8 #################################
# Load data ---------------------------------------------------------------
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/ManhattanPlots")
# data should be a csv file with the following column: locus, chrom, position
infile <- read.delim("batch_8_filteredMQ_filteredAS_aligned_loci.txt",header=TRUE)
head(infile)

fstfile <- read.delim("../../stacks_b8_wgenome_r05/orig_filter/batch_8_filteredMAF_filteredLoci30_filteredIndivids_filteredHWE_eastwest_fst_parsed.txt",header=TRUE)
head(fstfile)

# Load packages -----------------------------------------------------------
install.packages("gplots")
library(gplots)
install.packages("Hmisc")
library(Hmisc)
install.packages("plyr")
library(plyr)
install.packages("RColorBrewer")
library(RColorBrewer)
install.packages("dplyr")
library(dplyr)



# Join data frames --------------------------------------------------------
align_data <- left_join(infile,fstfile)
View(align_data)


# Testing with LG 13 -----------------------------------
## isolate only LG13, then join FST
LG13_df <- filter(align_data, LG == "LG13")
View(LG13_df)

## plot histogram
LG13_hist<-hist(LG13_df$Position, breaks=seq(0,255000000,150),
                xlab="Map Position (bp)",ylab="Number of PCod Loci", main = "Aligned PCod Loci Across ACod LG13")
length(LG13_df$Locus)
text(x = 2, y = 0, labels = "Total\nLoci:\n 188")


## Manhattan plot
install.packages("ggplot2")
library(ggplot2)

ggplot(LG13_df, aes(x=Position, y=Fst)) +
  geom_point()

ggplot(LG13_df, aes(x=Position, y=Fst)) +
  stat_bin_2d()


# Manhattan Plot of all Linkage Groups ----------------------------------------

ggplot(align_data, aes(x=Position, y=Fst)) +
  geom_line() +
  facet_wrap(~LG)

ggplot(align_data, aes(x=Position, y=Fst)) +
  geom_point() +
  facet_wrap(~LG)


# Manhattan Plots by Linkage Group ----------------------------------------
## LG11
LG11_df <- filter(align_data, LG == "LG11")
ggplot(LG11_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 11") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

ggplot(LG11_df, aes(x=Position, y=Fst)) +
  geom_point() +
  ggtitle("ACod Linkage Group 11") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))




## LG10
LG10_df <- filter(align_data, LG == "LG10")
ggplot(LG10_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 10") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG1
LG01_df <- filter(align_data, LG == "LG01")
ggplot(LG01_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 01") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG2
LG02_df <- filter(align_data, LG == "LG02")
ggplot(LG02_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 02") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))


## LG4
LG04_df <- filter(align_data, LG == "LG04")
ggplot(LG04_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 04") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))



## LG14
LG14_df <- filter(align_data, LG == "LG14")
ggplot(LG14_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 14") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))


## LG15
LG15_df <- filter(align_data, LG == "LG15")
ggplot(LG15_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 15") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG13
LG13_df <- filter(align_data, LG == "LG13")
ggplot(LG13_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 13") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG16
LG16_df <- filter(align_data, LG == "LG16")
ggplot(LG16_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 16") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG07
LG07_df <- filter(align_data, LG == "LG07")
ggplot(LG07_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 07") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))


## LG03
LG03_df <- filter(align_data, LG == "LG03")
ggplot(LG03_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 03") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))


## LG09
LG09_df <- filter(align_data, LG == "LG09")
ggplot(LG09_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 09") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))


## LG08
LG08_df <- filter(align_data, LG == "LG08")
ggplot(LG08_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 08") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG20
LG20_df <- filter(align_data, LG == "LG20")
ggplot(LG20_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 20") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG18
LG18_df <- filter(align_data, LG == "LG18")
ggplot(LG18_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 18") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG21
LG21_df <- filter(align_data, LG == "LG21")
ggplot(LG21_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 21") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG06
LG06_df <- filter(align_data, LG == "LG06")
ggplot(LG06_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 06") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG23
LG23_df <- filter(align_data, LG == "LG23")
ggplot(LG23_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 23") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG10
LG10_df <- filter(align_data, LG == "LG10")
ggplot(LG10_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 10") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG06
LG12_df <- filter(align_data, LG == "LG12")
ggplot(LG12_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 12") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG22
LG22_df <- filter(align_data, LG == "LG22")
ggplot(LG22_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 22") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))


## LG19
LG19_df <- filter(align_data, LG == "LG19")
ggplot(LG19_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 19") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG17
LG17_df <- filter(align_data, LG == "LG17")
ggplot(LG17_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 17") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))





############################  BATCH 4 #################################
# Load data ---------------------------------------------------------------
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/ManhattanPlots")
# data should be a csv file with the following column: locus, chrom, position
infile_b4 <- read.delim("batch_4_filteredMQ_filteredAS_aligned_loci.txt",header=TRUE)
head(infile)

fstfile_b4 <- read.delim("../../stacks_b4_wgenome/batch_4_MB_filteredMAF_filteredLoci50_filteredIndivids_filteredHWE_eastwest_fst_parsed.txt",header=TRUE)
head(fstfile)


# Join data frames --------------------------------------------------------
align_data_b4 <- left_join(infile_b4,fstfile_b4)
View(align_data_b4)


# Manhattan Plots by Linkage Group ----------------------------------------
## LG11
LG11_df <- filter(align_data, LG == "LG11")
ggplot(LG11_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 11") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

ggplot(LG11_df, aes(x=Position, y=Fst)) +
  geom_point() +
  ggtitle("ACod Linkage Group 11") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))




## LG10
LG10_df <- filter(align_data, LG == "LG10")
ggplot(LG10_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 10") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG1
LG01_df <- filter(align_data, LG == "LG01")
ggplot(LG01_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 01") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG2
LG02_df <- filter(align_data, LG == "LG02")
ggplot(LG02_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 02") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))


## LG4
LG04_df <- filter(align_data, LG == "LG04")
ggplot(LG04_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 04") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))



## LG14
LG14_df <- filter(align_data, LG == "LG14")
ggplot(LG14_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 14") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))


## LG15
LG15_df <- filter(align_data, LG == "LG15")
ggplot(LG15_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 15") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG13
LG13_df <- filter(align_data, LG == "LG13")
ggplot(LG13_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 13") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG16
LG16_df <- filter(align_data, LG == "LG16")
ggplot(LG16_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 16") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG07
LG07_df <- filter(align_data, LG == "LG07")
ggplot(LG07_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 07") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))


## LG03
LG03_df <- filter(align_data, LG == "LG03")
ggplot(LG03_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 03") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))


## LG09
LG09_df <- filter(align_data, LG == "LG09")
ggplot(LG09_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 09") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))


## LG08
LG08_df <- filter(align_data, LG == "LG08")
ggplot(LG08_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 08") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG20
LG20_df <- filter(align_data, LG == "LG20")
ggplot(LG20_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 20") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG18
LG18_df <- filter(align_data, LG == "LG18")
ggplot(LG18_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 18") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG21
LG21_df <- filter(align_data, LG == "LG21")
ggplot(LG21_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 21") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG06
LG06_df <- filter(align_data, LG == "LG06")
ggplot(LG06_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 06") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG23
LG23_df <- filter(align_data, LG == "LG23")
ggplot(LG23_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 23") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG10
LG10_df <- filter(align_data, LG == "LG10")
ggplot(LG10_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 10") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG06
LG12_df <- filter(align_data, LG == "LG12")
ggplot(LG12_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 12") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG22
LG22_df <- filter(align_data, LG == "LG22")
ggplot(LG22_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 22") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))


## LG19
LG19_df <- filter(align_data, LG == "LG19")
ggplot(LG19_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 19") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))

## LG17
LG17_df <- filter(align_data, LG == "LG17")
ggplot(LG17_df, aes(x=Position, y=Fst)) +
  geom_line() +
  ggtitle("ACod Linkage Group 17") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20))




# Manhattan Plot Across LGs -----------------------------------------------

posdata <- read.delim("batch_4_filteredMQ_filteredAS_aligned_loci_cpos.txt", header=TRUE)
View(posdata)

align_data_cpos <- left_join(posdata,fstfile_b4)
View(align_data_cpos)

x <- seq(0, 598315120, by=5000000)
y <- seq(0,1,by = 0.001)

ggplot(align_data_cpos, aes(x=C.Position, y=Fst)) + 
  stat_bin2d(breaks=list(x=x,y=y)) +
  ggtitle("ACod Linkage Groups Cumulative Position") + 
  theme(plot.title=element_text(hjust=0.5, size = 24)) + 
  theme(axis.title=element_text(size = 20)) +
  annotate("rect", xmin = 0, xmax = 28303952, ymin = 0, ymax = 1.0, alpha = 0.2) +
  annotate("rect", xmin = 52358358, xmax = 81809413, ymin = 0, ymax = 1.0, alpha = 0.2) +
  annotate("rect", xmin = 116614735, xmax = 140688790, ymin = 0, ymax = 1.0, alpha = 0.2) +
  annotate("rect", xmin = 166153410, xmax = 197386287, ymin = 0, ymax = 1.0, alpha = 0.2) +
  annotate("rect", xmin = 224183173, xmax = 249565487, ymin = 0, ymax = 1.0, alpha = 0.2) +
  annotate("rect", xmin = 274869793, xmax = 303812761, ymin = 0, ymax = 1.0, alpha = 0.2) +
  annotate("rect", xmin = 331110735, xmax = 356787470, ymin = 0, ymax = 1.0, alpha = 0.2) +
  annotate("rect", xmin = 386084402, xmax = 412682361, ymin = 0, ymax = 1.0, alpha = 0.2) +
  annotate("rect", xmin = 443775604, xmax = 462924811, ymin = 0, ymax = 1.0, alpha = 0.2) +
  annotate("rect", xmin = 485479066, xmax = 506655326, ymin = 0, ymax = 1.0, alpha = 0.2) +
  annotate("rect", xmin = 530804459, xmax = 553314763, ymin = 0, ymax = 1.0, alpha = 0.2) +
  annotate("rect", xmin = 575050466, xmax = 598315120, ymin = 0, ymax = 1.0, alpha = 0.2)
  
