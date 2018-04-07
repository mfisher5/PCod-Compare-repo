########### which linkage groups have the greatest number of islands? ############
#
# MF 4/6/2018
#
##################################################################################


# Install Packages --------------------------------------------------------
install.packages("dplyr")
install.packages("ggplot2")
library(readr)
library(ggplot2)
library(dplyr)



# Load Data -------------------------------------------------------------
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/SlidingWindow")

## east & west
mydata_islands <- read_delim("batch_8_eastwest_islands.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
head(mydata_islands)
colnames(mydata_islands) <- c("chromosome", "windows_per_island")

## west (note had to select out columns b/c i also had start & end positions in dataframe)
west_islands <- read_delim("batch_8_west_islands.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
head(west_islands)
colnames(west_islands) <- c("chromosome", "windows_per_island", "start", "end")
west_islands <- select(west_islands, chromosome, windows_per_island)
head(west_islands)


## east
east_islands <- read_delim("batch_8_east_islands.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
head(east_islands)
colnames(east_islands) <- c("chromosome", "windows_per_island", "start", "end")
east_islands <- select(east_islands, chromosome, windows_per_island)
head(east_islands)




# Calculate markers / islands per linkage group -------------------------------------
## count number of islands / markers per chromosome
mydata_islands_perlg <- mydata_islands %>%
  group_by(chromosome) %>%
  summarise(island_count = length(windows_per_island))
head(mydata_islands_perlg)



## west
west_islands_perlg <- west_islands %>%
  group_by(chromosome) %>%
  summarise(island_count = length(windows_per_island))
head(west_islands_perlg)




## east
east_islands_perlg <- east_islands %>%
  group_by(chromosome) %>%
  summarise(island_count = length(windows_per_island))
head(east_islands_perlg)



## all
islands_per_lg <- full_join(mydata_islands_perlg, west_islands_perlg, by = "chromosome")
head(islands_per_lg)
colnames(islands_per_lg) <- c("chromosome", "fulldata_count", "west_count")
head(islands_per_lg)
islands_per_lg <- full_join(islands_per_lg, east_islands_perlg, by="chromosome")
head(islands_per_lg)

# Plot islands per linkage group ------------------------------------------
library(reshape2)
islands_per_lg_melted <- melt(islands_per_lg)
islands_per_lg_melted

ggplot(islands_per_lg_melted, aes(x=chromosome, y=value, fill=variable)) +
  geom_col() +
  ylab("Number of Islands") +
  xlab("Linkage Group") +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=16), axis.title.y=element_text(size=16), axis.title.x = element_text(size=16), legend.text=element_text(size=16), legend.title = element_text(size=16)) +
  scale_fill_manual(values=c("grey48", "deepskyblue4", "mediumorchid2"), labels=c("All Data", "West", "East"), name="Analysis")
  





