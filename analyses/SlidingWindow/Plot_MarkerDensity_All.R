




# Install Packages --------------------------------------------------------

install.packages("ggplot2")
library(ggplot2)
library(dplyr)

# Load data ---------------------------------------------------------------
## set working directory
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/SlidingWindow")


all_marker_west <- read.table("West_SWA_Num_Name_Loci_Per_Window.txt", header=FALSE, sep="\t")
head(all_marker_west)
marker_west <- select(all_marker_west, c("V1", "V2", "V3"))
colnames(marker_west) <- c("chrom", "position", "num_markers")
head(marker_west)



all_marker_east <- read.table("East_Num_Name_Loci_Per_Window.txt", header=FALSE, sep="\t")
head(all_marker_east)
marker_east <- select(all_marker_east, c("V1", "V2", "V3"))
colnames(marker_east) <- c("chrom", "position", "num_markers")
head(marker_east)


all_marker <- read.table("EastWest_GlobalFst_Num_Name_Loci_Per_Window.txt", header=FALSE, sep="\t")
head(all_marker)
marker <- select(all_marker, c("V1", "V2", "V3"))
colnames(marker) <- c("chrom", "position", "num_markers")
head(marker)


marker <- mutate(marker, Analysis  = "All Data")
marker_east <- mutate(marker_east, Analysis = "East")
marker_west <- mutate(marker_west, Analysis = "West")

all_marker_data <- rbind(marker, marker_east, marker_west)



## create vector of number of markers per window

ggplot(all_marker_data, aes(x=num_markers)) +
  geom_histogram() + 
  stat_bin(binwidth=1) +
  facet_wrap(~Analysis) +
  xlab("Number of Markers per Window") +
  ylab("Number of Windows") +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=12))

qplot(markers_per_window, geom="histogram",
      binwidth = 1,
      main = "Markers per Sliding Window, West\nAcross All Linkage Groups",
      xlab = "# Markers", 
      ylab = "# Windows") +
  annotate("text", x = 25, y = 300, label="Mean:\n8.11")
