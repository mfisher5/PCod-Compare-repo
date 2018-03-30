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
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/stacks_pipeline_analyses/analyses_AK_LowGenotypeRate/batch8")

  
########################## PSTACKS

#Read in your datatable from a tab-delimited text file 
mydata <- read.delim("pstacks_loci_counts_b8_forR.txt")
#### HEADER SHOULD BE:(tab delimited)
#filename loci  population

View(mydata)

# Make the plot!!!


######### using facet_wrap to plot this
ggplot(mydata,aes(x=n_loci, fill = group)) + 
  geom_histogram(data = mydata, bins= 20) +
  facet_wrap(~population) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


######## produce the same plot, separated by population

mydata_kor <- filter(mydata, group == "Korea")
mydata_ak <- filter(mydata, group == "Alaska")
View(mydata_kor)

ggplot(mydata_kor,aes(x=n_loci, fill = group)) + 
  geom_histogram(data = mydata_kor, bins= 20, fill = "skyblue") +
  facet_wrap(~population) +
  xlab("Number of Loci Discovered") + 
  ylab("Number of Individuals") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12)) +
  xlim(0,27000)

ggplot(mydata_ak,aes(x=n_loci, fill = group)) + 
  geom_histogram(data = mydata_ak, bins= 20) +
  facet_wrap(~population) +
  xlab("Number of Loci Discovered") + 
  ylab("Number of Individuals") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=12)) +
  xlim(0,27000)





########################## SSTACKS

#Read in your datatable from a tab-delimited text file 
mydata_sstacks <- read.delim("sstacks_loci_counts_b8.txt")
#### HEADER SHOULD BE:(tab delimited)
#filename loci  population
View(mydata_sstacks)

# Make the plot!!!


######### using facet_wrap to plot this
ggplot(mydata_sstacks,aes(x=n_loci, fill = group)) + 
  geom_histogram(data = mydata_sstacks, bins= 20) +
  facet_wrap(~population) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



### 
mydata_kor_sstacks <- filter(mydata_sstacks, group == "Korea")
mydata_ak_sstacks <- filter(mydata_sstacks, group == "Alaska")
View(mydata_kor)

ggplot(mydata_kor_sstacks,aes(x=n_loci, fill = group)) + 
  geom_histogram(data = mydata_kor_sstacks, bins= 20, fill = "skyblue") +
  facet_wrap(~population) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(mydata_ak_sstacks,aes(x=n_loci, fill = group)) + 
  geom_histogram(data = mydata_ak_sstacks, bins= 20) +
  facet_wrap(~population) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlim(0,26000)




########################## SSTACKS DEPTH 10 +

#Read in your datatable from a tab-delimited text file 
mydata_d10 <- read.delim("sstacks_loci_counts_depth10_b8.txt")
#### HEADER SHOULD BE:(tab delimited)
#filename loci  population
View(mydata_d10)

# Make the plot!!!


######### using facet_wrap to plot this
ggplot(mydata_d10,aes(x=n_loci, fill = group)) + 
  geom_histogram(data = mydata_d10, bins= 20) +
  facet_wrap(~population) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



### 
mydata_kor_d10 <- filter(mydata_d10, group == "Korea")
mydata_ak_d10 <- filter(mydata_d10, group == "Alaska")

ggplot(mydata_kor_d10,aes(x=n_loci, fill = group)) + 
  geom_histogram(data = mydata_kor_d10, bins= 20, fill = "skyblue") +
  facet_wrap(~population) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(mydata_ak_d10,aes(x=n_loci, fill = group)) + 
  geom_histogram(data = mydata_ak_d10, bins= 20) +
  facet_wrap(~population) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlim(0,26000)



############################ GENEPOP, UNCORRECTED

#Read in your datatable from a tab-delimited text file 
mydata_pop <- read.delim("populations_loci_counts_b8.txt")
#### HEADER SHOULD BE:(tab delimited)
#filename loci  population
View(mydata_pop)

# Make the plot!!!


######### using facet_wrap to plot this
ggplot(mydata_pop,aes(x=n_loci, fill = group)) + 
  geom_histogram(data = mydata_pop, bins= 20) +
  facet_wrap(~population) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


############################ GENEPOP, CORRECTED

#Read in your datatable from a tab-delimited text file 
mydata_corr <- read.delim("populations_corrected_loci_counts_b8.txt")
#### HEADER SHOULD BE:(tab delimited)
#filename loci  population
View(mydata_corr)

# Make the plot!!!


######### using facet_wrap to plot this
ggplot(mydata_corr,aes(x=n_loci, fill = group)) + 
  geom_histogram(data = mydata_corr, bins= 20) +
  facet_wrap(~population) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

### 
mydata_kor_corr <- filter(mydata_corr, group == "Korea")
mydata_ak_corr <- filter(mydata_corr, group == "Alaska")

ggplot(mydata_kor_corr,aes(x=n_loci, fill = group)) + 
  geom_histogram(data = mydata_kor_corr, bins= 20, fill = "skyblue") +
  facet_wrap(~population) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(mydata_ak_corr,aes(x=n_loci, fill = group)) + 
  geom_histogram(data = mydata_ak_corr, bins= 20) +
  facet_wrap(~population) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlim(0,10000)


############################ GENEPOP, CORRECTED & FILTERED MAF


#Read in your datatable from a tab-delimited text file 
mydata_gen <- read.delim("genepop_loci_counts_b8_forR.txt")
#### HEADER SHOULD BE:(tab delimited)
#filename loci  population
View(mydata_gen)

# Make the plot!!!


######### using facet_wrap to plot this
ggplot(mydata_gen,aes(x=n_loci, fill = group)) + 
  geom_histogram(data = mydata_gen, bins= 20) +
  facet_wrap(~population) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

########################## SSTACKS DEPTH 5 +

#Read in your datatable from a tab-delimited text file 
mydata_d5 <- read.delim("sstacks_loci_counts_depth5_b8.txt")
#### HEADER SHOULD BE:(tab delimited)
#filename loci  population
View(mydata_d10)

# Make the plot!!!


######### using facet_wrap to plot this
ggplot(mydata_d5,aes(x=n_loci, fill = group)) + 
  geom_histogram(data = mydata_d5, bins= 20) +
  facet_wrap(~population) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



### 
mydata_kor_d5 <- filter(mydata_d5, group == "Korea")
mydata_ak_d5 <- filter(mydata_d5, group == "Alaska")

ggplot(mydata_kor_d5,aes(x=n_loci, fill = group)) + 
  geom_histogram(data = mydata_kor_d5, bins= 20, fill = "skyblue") +
  facet_wrap(~population) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(mydata_ak_d5,aes(x=n_loci, fill = group)) + 
  geom_histogram(data = mydata_ak_d5, bins= 20) +
  facet_wrap(~population) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlim(0,26000)












#### Use aggregate function to calculate
# summary statistics that will help you identify cut-off thresholds

# mean number of loci per pop
mymean <- aggregate(x = mydata$n_loci, by = list(mydata$population), mean)


### Use plyr package to estimate 1% quantile in the data per pop
my01quantile <- ddply(mydata, "population", summarise, pstacks_01_quantile = quantile(n_loci, .01))
mean(my01quantile$pstacks_01_quantile)

## Check up on a subset of individuals
trouble_samples <- mydata[mydata$n_loci<17390,]
good_samples <- mydata[mydata$n_loci>17390,]#Grab a subset of data with low numbers of loci present

# sort by number of loci in each sample (to help pick what individuals should go in the catalog)
sorted_samples <- good_samples[order(good_samples$n_loci),] 

write.table(trouble_samples, "samples_to_remove.txt", sep="\t")
write.table(good_samples, "samples_to_retain.txt", sep="\t")

write.table(sorted_samples, "samples_sorted_pstacks.txt", sep="\t")

# Extra notes:

  
# To find sweet colors in R: http://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3
# To find tutorials about ggplot: http://r-statistics.co/Top50-Ggplot2-Visualizations-MasterList-R-Code.html