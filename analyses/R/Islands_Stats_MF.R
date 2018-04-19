########### Distribution & T Test for Num/Width of Islands ############
#
# MF 4/19/2018
#
#######################################################################

library(ggplot2)
library(reshape2)
library(dplyr)

# Read in Data ------------------------------------------------------------
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/R")
num_data <- read.table("../../results/N_Islands_all.txt", sep="\t", header=TRUE)
width_data <- read.table("../../results/W_Islands_all.txt", sep="\t", header=TRUE)

head(num_data)
head(width_data)


# Plot --------------------------------------------------------------------
## melt data frame
num_melted <- melt(num_data)
colnames(num_melted) <- c("Data", "Number of Islands")
width_melted <- melt(width_data, na.rm=TRUE)
colnames(width_melted) <- c("Data", "Width of Islands")

## ggplot for number of islands
colors <- c("mediumorchid2", "deepskyblue4","grey48")

ggplot(num_melted, aes(x=Data, y=`Number of Islands`)) +
  geom_boxplot(fill=colors) +
  scale_x_discrete(labels=c("East", "West", "All Data")) +
  scale_y_continuous(breaks=seq(0,6,1)) +
  theme(axis.text = element_text(size=14), axis.title.y = element_text(size=16),
        axis.title.x = element_blank())


## ggplot for width of islands
colors <- c("mediumorchid2", "deepskyblue4","grey48")

ggplot(width_melted, aes(x=Data, y=`Width of Islands`)) +
  geom_boxplot(fill=colors) +
  scale_x_discrete(labels=c("East", "West", "All Data")) +
  scale_y_continuous(breaks=seq(0,12,1)) +
  theme(axis.text = element_text(size=14), axis.title.y = element_text(size=16),
        axis.title.x = element_blank())


# T-Test ------------------------------------------------------------------
width_data
hist(width_data$W.East); hist(width_data$W.West); hist(width_data$W.All)
hist(num_data$N.East); hist(num_data$N.West); hist(num_data$N.All)
## why can we do a t test if these are not normal distributions? The t-test is robust, especially if n1=n2 (it does) and the test is two-tailed (it is). Also, West and All variances are within 10% of each other

## WIDTH ##
t.test(x=width_data$W.East, y=width_data$W.West, alternative="two.sided", 
       mu=0, paired=FALSE, var.equal=FALSE, conf.level=0.95)
#Welch Two Sample t-test
#t = 0.48413, df = 51.217, p-value = 0.6304
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
# -0.8044368  1.3157868
#sample estimates:
#  mean of x mean of y 
#3.481481  3.225806 

t.test(x=width_data$W.East, y=width_data$W.All, alternative="two.sided", 
       mu=0, paired=FALSE, var.equal=FALSE, conf.level=0.95)
#Welch Two Sample t-test
#t = 2.0008, df = 37.151, p-value = 0.05276
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#   -0.01145095  1.83805027
#sample estimates:
#  mean of x mean of y 
#3.481481  2.568182 

t.test(x=width_data$W.West, y=width_data$W.All, alternative="two.sided", 
       mu=0, paired=FALSE, var.equal=TRUE, conf.level=0.95)
#Two Sample t-test
#t = 1.8475, df = 73, p-value = 0.06873
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.05179108  1.36704035
#sample estimates:
#  mean of x mean of y 
# 3.225806  2.568182 


## NUMBER ##

t.test(x=num_data$N.East, y=num_data$N.West, alternative="two.sided", 
       mu=0, paired=FALSE, var.equal=FALSE, conf.level=0.95)
#Welch Two Sample t-test
#t = -0.52035, df = 42.687, p-value = 0.6055
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.8480750  0.5002489
#sample estimates:
#  mean of x mean of y 
#1.173913  1.347826  

t.test(x=num_data$N.East, y=num_data$N.All, alternative="two.sided", 
       mu=0, paired=FALSE, var.equal=FALSE, conf.level=0.95)
#Welch Two Sample t-test
#t = -2.0217, df = 41.932, p-value = 0.04962
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#   -1.390108532 -0.001195816
#sample estimates:
#  mean of x mean of y 
#1.173913  1.869565

t.test(x=num_data$N.West, y=num_data$N.All, alternative="two.sided", 
       mu=0, paired=FALSE, var.equal=FALSE, conf.level=0.95)
#Welch Two Sample t-test
#t = -1.4045, df = 43.896, p-value = 0.1672
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
# -1.2704546  0.2269763
#sample estimates:
#  mean of x mean of y 
#1.347826  1.869565



# Kruskal-Wallace Test ---------------------------------------------------
kruskal.test(`Number of Islands` ~ Data, data=num_melted)
?kruskal.test

kruskal.test(`Width of Islands` ~ Data, data=width_melted)
