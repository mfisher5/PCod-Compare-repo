################# Manhattan Plot Code #################
#
# Modified from Charlie Waters' script to allow binning
#
# MF 3/25/2018
#
#######################################################

library(ggplot2)
library(dplyr)

# Load Data ---------------------------------------------------------------

## set working directory
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/ManhattanPlots")

## read in file with cumulative genome positions (created using excel)
align_data <- read.delim("../SlidingWindow/East/batch_8_SWA_input_east_sorted_cpos.txt",sep="\t",header=TRUE)
dim(align_data)
colnames(align_data)


align_data_west <- read.delim("../SlidingWindow/West/batch_8_SWA_input_west_sorted_cpos.txt",sep="\t",header=TRUE)
dim(align_data_west)
colnames(align_data_west)


## rename columns for R script function calls
colnames(align_data) <- c("Locus", "fst", "chromosome", "position")
View(align_data)

colnames(align_data_west) <- c("Locus", "fst", "chromosome", "position")



# Read in Outliers & Merge dataframes -------------------------------------
outliers_east <- read_csv("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/outliers/batch_8_final_filtered_aligned_EAST_outliers_2progs_forSWA.txt")

outlier_list <- c()
for(locus in align_data$Locus){
  if(locus %in% outliers_east$Locus){
    outlier_list <- c(outlier_list, "Outlier")
  } else {outlier_list <- c(outlier_list, "Neutral")}
}

plot_data_east <- cbind(align_data, outlier_list)
View(plot_data_east)




outliers_west <- read_csv("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/outliers/batch_8_final_filtered_aligned_WEST_outliers_2progs_forSWA.txt")

outlier_list_west <- c()
for(locus in align_data_west$Locus){
  if(locus %in% outliers_west$Locus){
    outlier_list_west <- c(outlier_list_west, "Outlier")
  } else {outlier_list_west <- c(outlier_list_west, "Neutral")}
}

plot_data_west <- cbind(align_data_west, outlier_list_west)
View(plot_data_west)




# Create bins -------------------------------------------------------------
## create bin breaks using dplyr

plot_data_east <- mutate(plot_data_east, bin = ntile(position, 200))
plot_data_west <- mutate(plot_data_west, bin = ntile(position, 200))

## create actual bin values
newbins = seq(0, max(plot_data_east$position), length.out = 200)
head(newbins)

newbins_west = seq(0, max(plot_data_west$position), length.out = 200)
head(newbins_west)

## use dplyr bins to create vector of actual bin values
newbinlist <- c()

for(i in seq(1,200)){
  bindata <- filter(plot_data_east, bin == i)
  newbin <- max(bindata$position)
  newbinlist <- c(newbinlist, rep(newbin, times = length(bindata$position)))
}


newbinlist_west <- c()

for(i in seq(1,200)){
  bindata <- filter(plot_data_west, bin == i)
  newbin <- max(bindata$position)
  newbinlist_west <- c(newbinlist_west, rep(newbin, times = length(bindata$position)))
}


## check the number of bin values created is = length of data frame
length(newbinlist)
length(newbinlist_west)


## remove dplyr bins and add in new bins to data frame
plot_data_east_binned <- plot_data_east %>%
  select(c("Locus", "fst", "chromosome","position", "outlier_list")) %>%
  mutate(bin = newbinlist)
View(plot_data_east_binned)

length(plot_data_east_binned$Locus)



plot_data_west_binned <- plot_data_west %>%
  select(c("Locus", "fst", "chromosome","position", "outlier_list_west")) %>%
  mutate(bin = newbinlist_west)
View(plot_data_west_binned)

length(plot_data_west_binned$Locus)




# Extract Outliers; Make 1 Data Frame -------------------------------------
plot_east_outliers <- plot_data_east_binned %>%
  filter(outlier_list == "Outlier") %>%
  mutate(Population = "East")
plot_west_outliers <- plot_data_west_binned %>%
  filter(outlier_list_west == "Outlier") %>%
  mutate(Population = "West")
head(plot_west_outliers)

colnames(plot_west_outliers) <- c("Locus","fst","chromosome","position","outlier_list","bin","Population")

plot_outliers_data <- rbind(plot_east_outliers, plot_west_outliers)


# Plot binned Manhattan Plot ----------------------------------------------
ggplot(plot_outliers_data, aes(x = bin, y = fst)) +
  geom_point(aes(col = factor(Population), pch = factor(Population)), size = 4) +
  ggtitle("Outlier Loci\nPer Locus Fst") + 
  xlab("Linkage Group") +
  ylab(expression("F"[st])) +
  theme(plot.title=element_text(hjust=0.5, size = 16)) + 
  theme(axis.title=element_text(size = 16)) +
  theme(axis.text.x = element_blank(), axis.text.y=element_text(size=14)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.background = element_blank()) +
  scale_color_manual(values = c("mediumorchid2", "deepskyblue4")) +
  scale_y_continuous(limits=(c(-0.02,0.5)), expand=c(0,0)) +
  scale_x_continuous(limits=c(0,598315121), expand=c(0,0)) +
  annotate("rect", xmin = 0, xmax = 28303952, ymin = 0, ymax = 0.5, alpha = 0.2) +
  annotate("rect", xmin = 52358358, xmax = 81809413, ymin = 0, ymax = 0.5, alpha = 0.2) +
  annotate("rect", xmin = 116614735, xmax = 140688790, ymin = 0, ymax = 0.5, alpha = 0.2) +
  annotate("rect", xmin = 166153410, xmax = 197386287, ymin = 0, ymax = 0.5, alpha = 0.2) +
  annotate("rect", xmin = 224183173, xmax = 249565487, ymin = 0, ymax = 0.5, alpha = 0.2) +
  annotate("rect", xmin = 274869793, xmax = 303812761, ymin = 0, ymax = 0.5, alpha = 0.2) +
  annotate("rect", xmin = 331110735, xmax = 356787470, ymin = 0, ymax = 0.5, alpha = 0.2) +
  annotate("rect", xmin = 386084402, xmax = 412682361, ymin = 0, ymax = 0.5, alpha = 0.2) +
  annotate("rect", xmin = 443775604, xmax = 462924811, ymin = 0, ymax = 0.5, alpha = 0.2) +
  annotate("rect", xmin = 485479066, xmax = 506655326, ymin = 0, ymax = 0.5, alpha = 0.2) +
  annotate("rect", xmin = 530804459, xmax = 553314763, ymin = 0, ymax = 0.5, alpha = 0.2) +
  annotate("rect", xmin = 575050466, xmax = 598315120, ymin = 0, ymax = 0.5, alpha = 0.2)




