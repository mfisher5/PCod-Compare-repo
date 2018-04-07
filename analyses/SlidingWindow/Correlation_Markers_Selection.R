
# Install Packages --------------------------------------------------------
install.packages("dplyr")
install.packages("ggplot2")
library(readr)
library(ggplot2)
library(dplyr)



# Load Data -------------------------------------------------------------
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/SlidingWindow")

## east & west
all_marker_data <- read.table("EastWest_GlobalFst_Num_Name_Loci_Per_Window.txt", header=FALSE, sep="\t")
head(all_marker_data)
mydata_markers <- select(all_marker_data, c("V1", "V2", "V3"))
colnames(mydata_markers) <- c("chromosome", "position", "num_markers")

eastwest <- read_delim("EastvWest/batch_8_SWA_eastwest_output_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
head(eastwest)

eastwest_selection <- eastwest %>%
  mutate(positive = ifelse(`Fst/Fct` > upper_95, "1", "0")) %>%
  mutate(negative = ifelse(`Fst/Fct` < lower_95, "1", "0")) %>%
  mutate(selection = ifelse(positive == 1 | negative == 1, "Yes", "No"))
head(eastwest_selection)


mydata_joined <- full_join(eastwest_selection, mydata_markers, by="position")
head(mydata_joined)

mydata_joined <- select(mydata_joined, c(chromosome.x, position, selection, num_markers))
mydata_joined_noNA <- filter(mydata_joined, is.na(selection) == FALSE)
head(mydata_joined_noNA)

ggplot(mydata_joined_noNA, aes(x=selection, y=num_markers, fill = selection)) +
  geom_boxplot() +
  xlab("Does Window exceed 95% Confidence Interval?") +
  ylab("Number of Markers per Window") +
  theme(axis.title=element_text(size=16), axis.text.y = element_text(size=16))



## west

all_marker_data <- read.table("West_SWA_Num_Name_Loci_Per_Window.txt", header=FALSE, sep="\t")
head(all_marker_data)
west_markers <- select(all_marker_data, c("V1", "V2", "V3"))
colnames(west_markers) <- c("chromosome", "position", "num_markers")

west <- read_delim("West/batch_8_final_filtered_west_2reg_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
head(west)

west_selection <- west %>%
  mutate(positive = ifelse(`Fst/Fct` > upper_95, "1", "0")) %>%
  mutate(negative = ifelse(`Fst/Fct` < lower_95, "1", "0")) %>%
  mutate(selection = ifelse(positive == 1 | negative == 1, "Yes", "No"))
head(west_selection)


west_joined <- full_join(west_selection, west_markers, by="position")
head(west_joined)

west_joined <- select(west_joined, c(chromosome.x, position, selection, num_markers))
head(west_joined)
west_joined_noNA <- filter(west_joined, is.na(selection) == FALSE)
head(west_joined_noNA)

ggplot(west_joined_noNA, aes(x=selection, y=num_markers, fill = selection)) +
  geom_boxplot() +
  xlab("Does Window exceed 95% Confidence Interval?") +
  ylab("Number of Markers per Window") +
  theme(axis.title=element_text(size=16), axis.text.y = element_text(size=16))





## east

all_marker_data <- read.table("East_Num_Name_Loci_Per_Window.txt", header=FALSE, sep="\t")
head(all_marker_data)
east_markers <- select(all_marker_data, c("V1", "V2", "V3"))
colnames(east_markers) <- c("chromosome", "position", "num_markers")

east <- read_delim("East/batch_8_final_filtered_east_globalFst_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
head(east)

east_selection <- east %>%
  mutate(positive = ifelse(`Fst/Fct` > upper_95, "1", "0")) %>%
  mutate(negative = ifelse(`Fst/Fct` < lower_95, "1", "0")) %>%
  mutate(selection = ifelse(positive == 1 | negative == 1, "Yes", "No"))
head(east_selection)


east_joined <- full_join(east_selection, east_markers, by="position")
head(east_joined)

east_joined <- select(east_joined, c(chromosome.x, position, selection, num_markers))
head(east_joined)
east_joined_noNA <- filter(east_joined, is.na(selection) == FALSE)
head(east_joined_noNA)

ggplot(east_joined_noNA, aes(x=selection, y=num_markers, fill = selection)) +
  geom_boxplot() +
  xlab("Does Window exceed 95% Confidence Interval?") +
  ylab("Number of Markers per Window") +
  theme(axis.title=element_text(size=16), axis.text.y = element_text(size=16))








