############### Plot PCA: neutral and outlier loci #################
#
# MF 4/4/2018
# PCod-Compare
#
####################################################################


# Load Packages -----------------------------------------------------------
install.packages("adegenet")
install.packages("hierfstat")
install.packages("gplots")
library(adegenet)
library(hierfstat)
library (gplots)


# Read in Data as a genepop file ------------------------------------------------------------
# File can be delimited by tabs or spaces but there must abe a comma after each individual. 
# Specify how many characters code each allele with ncode. 

setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/PCA")

my_data_outliers <-read.genepop("batch_8_final_filtered_aligned_genepop_eastwest_outliers.gen")

my_data_2outliers <- read.genepop("batch_8_final_filtered_aligned_genepop_eastwest_outliers_2progs.gen")

my_data_neutral <- read.genepop("batch_8_final_filtered_aligned_genepop_eastwest_neutral.gen")


# To retreive useful data summaries
(summary(my_data_outliers))
my_data_outliers$pop

(summary(my_data_2outliers))

(summary(my_data_neutral))

# To replace missing data information with the mean
X <- scaleGen(my_data_outliers, NA.method="mean")
X.2o <- scaleGen(my_data_2outliers, NA.method="mean")
X.n <- scaleGen(my_data_neutral, NA.method="mean")




# Run PCA -----------------------------------------------------------------
# To conduct the PCA. IF YOU DO NOT KNOW HOW MANY AXES TO RETAIN
## for all outliers
pca_outliers <- dudi.pca(X,cent=FALSE,scale=FALSE)
barplot(pca_outliers$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
summary(pca_outliers)

## for all outliers in 2 programs or more
pca_2outliers <- dudi.pca(X.2o,cent=FALSE,scale=FALSE)
barplot(pca_2outliers$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

## for all neutral loci
pca_neutral <- dudi.pca(X.n,cent=FALSE,scale=FALSE)
barplot(pca_neutral$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))



# To conduct the PCA. IF YOU KNOW HOW MANY AXES TO RETAIN
pca_outliers <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=2) #nf = number to retain
barplot(pca_outliers$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

pca_2outliers <- dudi.pca(X.2o,cent=FALSE,scale=FALSE,scannf=FALSE,nf=2) #nf = number to retain
barplot(pca_2outliers$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

pca_neutral <- dudi.pca(X.n,cent=FALSE,scale=FALSE,scannf=FALSE,nf=2) #nf = number to retain
barplot(pca_outliers$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))




# Plot Outlier PCAs ----------------------------------------------------------------

#Color of points (order in genepop)
col <- c("deepskyblue4","deepskyblue4","deepskyblue4","green4","deepskyblue4","green4","darkorchid2","darkorchid2","darkorchid2","darkorchid2","darkorchid2","darkorchid2")
col_leg <- c("deepskyblue4","green4","darkorchid2")

#Point shapes for legend (here is by region)
points_leg <- c(16,16,15)

#to graph with lines between samples
s.class(pca_outliers$li, fac=pop(my_data_outliers), 
        col=col, #color of points. will retain lines between points
        clabel=0, #remove population labels
        cellipse=0, #remove ellipses; to add back in, make >=1
        cpoint=1,
        grid=FALSE, #otherwise will have light gray grid markers
        pch=c(16,16,16,16,16,16,15,15,15,15,15,15)[as.numeric(pop(my_data_outliers))], #change point shapes
        axesell=TRUE)
# add legend
legend (x=10,y=25, legend = c("West - Southern Coast", "West - Western Coast", "East"), col = col_leg, border = FALSE, bty = "n", cex = 0.9, pt.cex=1.5, title = "Population",pch=points_leg)
# add eigenvalues plot as inset
add.scatter.eig(pca_outliers$eig[1:50],posi="bottomright", 3,2,1,ratio=.3)
# percent variation explained
eig.perc <- 100*pca_outliers$eig/sum(pca_outliers$eig)
head(eig.perc)




#to graph with lines between samples
s.class(pca_2outliers$li, fac=pop(my_data_2outliers), 
        col=col, #color of points. will retain lines between points
        clabel=0, #remove population labels
        cellipse=0, #remove ellipses; to add back in, make >=1
        cpoint=1,
        grid=FALSE, #otherwise will have light gray grid markers
        pch=c(16,16,16,16,16,16,15,15,15,15,15,15)[as.numeric(pop(my_data_2outliers))], #change point shapes
        axesell=TRUE)
# add legend
legend (x=-16,y=14, legend = c("West-Southern Coast", "West-Western Coast", "East"), col = col_leg, border = FALSE, bty = "n", cex = 0.9, pt.cex=1.5, title = "Population",pch=points_leg)
# add eigenvalues plot as inset
add.scatter.eig(pca_2outliers$eig[1:50],posi="bottom", 3,2,1,ratio=.3)
# percent variation explained
eig.perc <- 100*pca_2outliers$eig/sum(pca_2outliers$eig)
head(eig.perc)




# Plot Neutral PCAs ----------------------------------------------------------------

#Color of points (order in genepop)
col <- c("deepskyblue4","deepskyblue4","deepskyblue4","green4","deepskyblue4","green4","darkorchid2","darkorchid2","darkorchid2","darkorchid2","darkorchid2","darkorchid2")
col_leg <- c("deepskyblue4","green4","darkorchid2")

#Point shapes for legend (here is by region)
points_leg <- c(16,16,15)

#to graph with lines between samples
s.class(pca_neutral$li, fac=pop(my_data_neutral), 
        col=col, #color of points. will retain lines between points
        clabel=0, #remove population labels
        cellipse=0, #remove ellipses; to add back in, make >=1
        cpoint=1,
        grid=FALSE, #otherwise will have light gray grid markers
        pch=c(16,16,16,16,16,16,15,15,15,15,15,15)[as.numeric(pop(my_data_neutral))], #change point shapes
        axesell=TRUE)
# add legend
legend (x=-2.8,y=25, legend = c("West - Southern Coast", "West - Western Coast", "East"), col = col_leg, border = FALSE, bty = "n", cex = 0.9, pt.cex=1.5, title = "Population",pch=points_leg)
# add eigenvalues plot as inset
add.scatter.eig(pca_neutral$eig[1:50],posi="topright", 3,2,1,ratio=.3)
# percent variation explained
eig.perc <- 100*pca_neutral$eig/sum(pca_neutral$eig)
head(eig.perc)



