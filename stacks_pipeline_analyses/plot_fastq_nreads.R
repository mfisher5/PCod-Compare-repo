# Script by Eleni Petrou. December 27, 2017


# Load necessary libraries
install.packages("ggalt")
library(ggplot2)
library(ggalt)
library(RColorBrewer)
library(plyr)
library(dplyr)

# Specify that you don't want scientific notation
options(scipen = 999) #Positive values bias towards fixed and negative towards scientific notation

# Set the bw theme (layout) for ggplot.
theme_set(theme_bw()) 


# Set working directory 
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/stacks_pipeline_analyses")



########################## PSTACKS

#Read in your datatable from a tab-delimited text file 
mydata <- read.delim("fastq_readcounts_forR.txt")
#### HEADER SHOULD BE:(tab delimited)
#filename loci  population

head(mydata)

# Make the plot!!!


######### using facet_wrap to plot this
ggplot(mydata,aes(x=n_reads, fill = group)) + 
  geom_histogram(data = mydata, bins= 20) +
  facet_wrap(~population) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


######## produce the same plot, separated by population

mydata_kor <- filter(mydata, group == "Korea")
mydata_kor <- mutate(mydata_kor, th_reads = n_reads/100000)
head(mydata_kor)
mydata_ak <- mydata %>%
  filter(group == "Alaska") %>%
  mutate(th_reads = n_reads/100000)
head(mydata_ak)



ggplot(mydata_kor,aes(x=th_reads, fill = group)) + 
  geom_histogram(data = mydata_kor, bins= 20, fill = "skyblue") +
  facet_wrap(~population) +
  xlab("Number of Raw Reads (x 100,000)") + 
  ylab("Number of Individuals") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12)) +
  xlim(0,200)

ggplot(mydata_ak,aes(x=th_reads, fill = group)) + 
  geom_histogram(data = mydata_ak, bins= 20) +
  facet_wrap(~population) +
  xlab("Number of Raw Reads (x 100,000)") + 
  ylab("Number of Individuals") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12)) +
  xlim(0,200)


