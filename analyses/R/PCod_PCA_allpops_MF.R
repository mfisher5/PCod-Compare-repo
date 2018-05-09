############### Plot PCA: neutral and outlier loci #################
#
# MF 4/4/2018
# PCod-Compare
# Most code borrowed from Eleni Petrou
# Edited 5/7/2018 to display all three Korean populations of PCod
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
mydata <- read.genepop("../../stacks_b8_wgenome_r05/batch_8_final_filtered_aligned_genepop.gen")

# To retreive useful data summaries
(summary(mydata))


# To replace missing data information with the mean
X <- scaleGen(mydata, NA.method="mean")




# Run PCA -----------------------------------------------------------------
# To conduct the PCA. IF YOU DO NOT KNOW HOW MANY AXES TO RETAIN\
pca_mydata <- dudi.pca(X,cent=FALSE,scale=FALSE)
barplot(pca_mydata$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
summary(pca_mydata)



# To conduct the PCA. IF YOU KNOW HOW MANY AXES TO RETAIN
pca_mydata <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=2) #nf = number to retain
barplot(pca_mydata$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))



# Plot Outlier PCAs ----------------------------------------------------------------

#Color of points (by region; to find order in genepop, use levels(pop(mydata)))
col <- c("deepskyblue4","deepskyblue4","deepskyblue4","aquamarine3","chartreuse4","deepskyblue4","deepskyblue4","aquamarine3","deepskyblue4","darkorchid2","darkorchid2","darkorchid2","darkorchid2","darkorchid2","darkorchid2")
col_leg <- c("deepskyblue4","aquamarine3","chartreuse4","darkorchid2")

#Point shapes for legend (here is by region)
points_leg <- c(16,16,16,15)

#to graph with lines between samples
s.class(pca_mydata$li, fac=pop(mydata), 
        col=col, #color of points. will retain lines between points
        clabel=0, #remove population labels
        cellipse=0, #remove ellipses; to add back in, make >=1
        cpoint=1,
        grid=FALSE, #otherwise will have light gray grid markers
        pch=c(16,16,16,16,16,16,16,16,16,15,15,15,15,15,15,15,15)[as.numeric(pop(mydata))], #change point shapes
        axesell=TRUE)
# add legend
legend (x=10,y=40, legend = c("West - Southern Coast", "West - Western Coast", "West - Eastern Coast","East"), col = col_leg, border = FALSE, bty = "n", cex = 1.2, pt.cex=1.5, title = "Population",pch=points_leg)
# add eigenvalues plot as inset
add.scatter.eig(pca_mydata$eig[1:50],posi="bottomright", 3,2,1,ratio=.3)
# percent variation explained
eig.perc <- 100*pca_mydata$eig/sum(pca_mydata$eig)
head(eig.perc)

