########### Correlate # of islands to # of markers ############
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


mydata_markers <- read_delim("EastvWest/batch_8_SWA_input_eastwest_globalFst_filtered_sorted.txt","\t", escape_double = FALSE, trim_ws = TRUE)
head(mydata_markers)


## west (note had to select out columns b/c i also had start & end positions in dataframe)
west_islands <- read_delim("batch_8_west_islands.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
head(west_islands)
colnames(west_islands) <- c("chromosome", "windows_per_island", "start", "end")
west_islands <- select(west_islands, chromosome, windows_per_island)
head(west_islands)

west_markers <- read_delim("West/batch_8_SWA_input_west_sorted.txt","\t", escape_double = FALSE, trim_ws = TRUE)
head(west_markers)




## east
east_islands <- read_delim("batch_8_east_islands.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
head(east_islands)
colnames(east_islands) <- c("chromosome", "windows_per_island", "start", "end")
east_islands <- select(east_islands, chromosome, windows_per_island)
head(east_islands)


east_markers <- read_delim("East/batch_8_SWA_input_east_sorted.txt","\t", escape_double = FALSE, trim_ws = TRUE)
head(east_markers)





# Calculate islands per linkage group -------------------------------------
## count number of islands / markers per chromosome
mydata_islands_perlg <- mydata_islands %>%
  group_by(chromosome) %>%
  summarise(island_count = length(windows_per_island)) %>%
  rbind(c("LG03", 0)) %>%
  rbind(c("LG05", 0)) %>%
  rbind(c("LG12", 0)) %>%
  rbind(c("LG22", 0))
tail(mydata_islands_perlg)




## west
west_islands_perlg <- west_islands %>%
  group_by(chromosome) %>%
  summarise(island_count = length(windows_per_island)) %>%
  rbind(c("LG05", 0)) %>%
  rbind(c("LG06", 0)) %>%
  rbind(c("LG07", 0)) %>%
  rbind(c("LG18", 0)) %>%
  rbind(c("LG23", 0))
tail(west_islands_perlg)

unique(west_islands_perlg$chromosome)




## east
east_islands_perlg <- east_islands %>%
  group_by(chromosome) %>%
  summarise(island_count = length(windows_per_island)) %>%
  rbind(c("LG01", 0)) %>%
  rbind(c("LG08", 0)) %>%
  rbind(c("LG11", 0)) %>%
  rbind(c("LG18", 0)) %>%
  rbind(c("LG20", 0))
tail(east_islands_perlg)

unique(east_islands_perlg$chromosome)



# Calculate # markers per LG ----------------------------------------------

mydata_markers_perlg <- mydata_markers %>%
  group_by(chromosome) %>%
  summarise(marker_count = length(Locus))
head(mydata_markers_perlg)

west_markers_perlg <- west_markers %>%
  group_by(chromosome) %>%
  summarise(marker_count = length(Locus))
head(west_markers_perlg)


east_markers_perlg <- east_markers %>%
  group_by(chromosome) %>%
  summarise(marker_count = length(Locus))
head(east_markers_perlg)





# Calculate width per island, join ----------------------------------------

mydata_width_perlg <- mydata_islands %>%
  group_by(chromosome) %>%
  summarise(island_width = mean(windows_per_island))
head(mydata_width_perlg)


## join data frames by loci (left join b/c not all chromosomes have islands)
mydata_islands_perlg$chromosome <- as.character(mydata_islands_perlg$chromosome)
mydata_width_perlg$chromosome <- as.character(mydata_width_perlg$chromosome)
mydata_joined <- left_join(mydata_islands_perlg, mydata_markers_perlg, by="chromosome")
mydata_joined <- full_join(mydata_joined, mydata_width_perlg, by="chromosome")
mydata_joined


## west
west_width_perlg <- west_islands %>%
  group_by(chromosome) %>%
  summarise(island_width = mean(windows_per_island))
head(west_width_perlg)

west_islands_perlg$chromosome <- as.character(west_islands_perlg$chromosome)
west_width_perlg$chromosome <- as.character(west_width_perlg$chromosome)
west_joined <- left_join(west_islands_perlg, west_markers_perlg, by="chromosome")
west_joined <- full_join(west_joined, west_width_perlg, by="chromosome")
west_joined

## east
east_width_perlg <- east_islands %>%
  group_by(chromosome) %>%
  summarise(island_width = mean(windows_per_island))
head(east_width_perlg)

east_islands_perlg$chromosome <- as.character(east_islands_perlg$chromosome)
east_width_perlg$chromosome <- as.character(east_width_perlg$chromosome)
east_joined <- left_join(east_islands_perlg, east_markers_perlg, by="chromosome")
east_joined <- full_join(east_joined, east_width_perlg, by="chromosome")
east_joined




# Plot Number Islands v. Number Markers -----------------------------------

fit <- lm(island_count ~ marker_count, data=mydata_joined)
summary(fit)

ggplot(mydata_joined, aes(x=marker_count,y=island_count)) +
  geom_point(size=2) +
  xlab("Number of Markers") +
  ylab("Number of Islands") +
  ggtitle("All Data") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20)) +
  annotate("text", label="R-squared:0.127; p = 0.053", x=150,y=5, color="lightsteelblue4", size=6)


## west
fit <- lm(island_count ~ marker_count, data=west_joined)
summary(fit)

ggplot(west_joined, aes(x=marker_count,y=island_count)) +
  geom_point(size=2) +
  xlab("Number of Markers") +
  ylab("Number of Islands") +
  ggtitle("West") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20)) +
  annotate("text", label="R-squared:0.00947; p = 0.2839", x=110,y=5, color="lightsteelblue4", size=6)

## east
fit <- lm(island_count ~ marker_count, data=east_joined)
summary(fit)

ggplot(east_joined, aes(x=marker_count,y=island_count)) +
  geom_point(size=2) +
  xlab("Number of Markers") +
  ylab("Number of Islands") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=20)) +
  annotate("text", label="R-squared:-0.02328, p = 0.488", x=110,y=5, color="lightsteelblue4", size = 6)



## all

fit <- lm(island_count ~ marker_count, data=all_marker_island_data)
summary(fit)
ggplot(all_marker_island_data, aes(x=marker_count,y=island_count)) +
  geom_point(size=2) +
  xlab("Number of Markers") +
  ylab("Number of Islands") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=20)) +
  geom_smooth(method="lm", se=F, linetype = 2, color="lightsteelblue4", size=0.6) +
  annotate("text", label="R-squared:0.2364; p value = 9.97e-5", x=150,y=4.5, color="black", size = 6)




