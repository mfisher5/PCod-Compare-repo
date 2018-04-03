################# Manhattan Plot Code #################
#
# Modified from Charlie Waters' script to allow binning
#
# MF 3/25/2018
#
#######################################################

library(ggplot2)

# Load Data ---------------------------------------------------------------

## set working directory
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/ManhattanPlots")

## read in file with cumulative genome positions (created using excel)
align_data <- read.delim("../SlidingWindow/West/batch_8_SWA_input_west_sorted_cpos.txt",sep="\t",header=TRUE)
dim(align_data)
colnames(align_data)

## rename columns for R script function calls
colnames(align_data) <- c("Locus", "fst", "chromosome", "position")
View(align_data)




# Create bins -------------------------------------------------------------
## create bin breaks using dplyr
library(dplyr)

align_data <- mutate(align_data, bin = ntile(position, 200))

## create actual bin values
newbins = seq(0, max(align_data$position), length.out = 200)
head(newbins)

## use dplyr bins to create vector of actual bin values
newbinlist <- c()

for(i in seq(1,200)){
  bindata <- filter(align_data, bin == i)
  newbin <- max(bindata$position)
  newbinlist <- c(newbinlist, rep(newbin, times = length(bindata$position)))
}

## check the number of bin values created is = length of data frame
length(newbinlist)


## remove dplyr bins and add in new bins to data frame
align_data_binned <- align_data %>%
  select(c("Locus", "fst", "chromosome","position")) %>%
  mutate(bin = newbinlist)
View(align_data_binned)




# Plot binned Manhattan Plot ----------------------------------------------
ggplot(align_data_binned, aes(x=bin, y=fst)) + 
  geom_point(size=1) +
  ggtitle("Per Locus Fst, West Pacific cod") + 
  xlab("Atlantic cod Linkage Group") +
  ylab(expression("F"[st])) +
  theme(plot.title=element_text(hjust=0.5, size = 20)) + 
  theme(axis.title=element_text(size = 16)) +
  theme(axis.text.x = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.background = element_blank()) +
  scale_y_continuous(limits=(c(-0.02,1)), expand=c(0,0)) +
  scale_x_continuous(limits=c(0,598315121), expand=c(0,0)) +
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


598315120/200




# Highlight Outliers ------------------------------------------------------
## add column to data frame which indicates outlier loci
outliers <- read.table("batch_8_final_filtered_aligned_WEST_outliers.txt", header=TRUE)
head(outliers)

align_data_binned_outliers <- left_join(align_data_binned, outliers, by="Locus")
align_data_binned_outliers$Program[is.na(align_data_binned_outliers$Program)] <- 0
View(align_data_binned_outliers)

## plot binned Manhattan plot with outliers highlighted

ggplot(align_data_binned_outliers, aes(x=bin, y=fst)) + 
  geom_point(aes(colour=factor(Program), size=factor(Program))) +
  scale_colour_manual(values=c("black", "deeppink2", "deeppink4"), labels=c("None", "OutFLANK", "OutFLANK+Bayescan"), name = "Outliers") +
  scale_size_manual(values=c(1,2,2), labels=c("None", "OutFLANK", "OutFLANK+Bayescan"), name = "Outliers") +
  ggtitle("Per Locus Fst, West Pacific cod") + 
  xlab("Atlantic cod Linkage Group") +
  ylab(expression("F"[st])) +
  theme(plot.title=element_text(hjust=0.5, size = 20)) + 
  theme(axis.title=element_text(size = 16)) +
  theme(axis.text.x = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.background = element_blank()) +
  scale_y_continuous(limits=(c(-0.02,1)), expand=c(0,0)) +
  scale_x_continuous(limits=c(0,598315121), expand=c(0,0)) +
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

align_data_binned_outliers_only <- filter(align_data_binned_outliers, Program == 1 | Program == 3)
View(align_data_binned_outliers_only)
