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

my_data_outliers <-read.genepop("batch_8_final_filtered_aligned_genepop_east_coastal_outliers.gen")

my_data_2outliers <- read.genepop("batch_8_final_filtered_aligned_genepop_east_coastal_outliers2p.gen")

my_data_neutral <- read.genepop("batch_8_final_filtered_aligned_genepop_east_coastal_neutral.gen")


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
barplot(pca_neutral$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))




# Plot Outlier PCAs ----------------------------------------------------------------

#Color of points (order in genepop)
col <- c("tomato1","#bd0026","gold", "#feb24c", "sienna", "orangered2")
names_leg <- c("Wash. Coast", "Hecate Strait", "PW Sound", "Kodiak", "Unimak Pass", "Adak")
col_leg <- c("gold", "#feb24c", "sienna", "tomato1", "orangered2", "#bd0026")

#to graph with lines between samples
s.class(pca_outliers$li, fac=pop(my_data_outliers), 
        col=col, #color of points. will retain lines between points
        clabel=0, #remove population labels
        cellipse=0, #remove ellipses; to add back in, make >=1
        cpoint=1,
        grid=FALSE, #otherwise will have light gray grid markers
        pch=19, #change point shapes
        axesell=TRUE)
# add legend
legend(x=-13,y=7.5, legend = names_leg, col = col_leg, border = FALSE, bty = "n", cex = 0.9, pt.cex=1.5, title = "Sampling Site",pch=19)
# add eigenvalues plot as inset
add.scatter.eig(pca_outliers$eig[1:50],posi="bottomright", 3,2,1,ratio=.3)

# check legend colors
legend("center", legend = unique(pop(my_data_outliers)), col = col, pt.cex=1.5, pch = 19)
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
legend(x=-6.5,y=5, legend = names_leg, col = col_leg, border = FALSE, bty = "n", cex = 0.9, pt.cex=1.5, title = "Sampling Site",pch=19)
# add eigenvalues plot as inset
add.scatter.eig(pca_2outliers$eig[1:25],posi="bottomright", 3,2,1,ratio=.3)
# percent variation explained
eig.perc <- 100*pca_2outliers$eig/sum(pca_2outliers$eig)
head(eig.perc)



# Plot Neutral PCAs ----------------------------------------------------------------

#Color of points (order in genepop)
col <- c("tomato1","#bd0026","gold", "#feb24c", "sienna", "orangered2")
names_leg <- c("Wash. Coast", "Hecate Strait", "PW Sound", "Kodiak", "Unimak Pass", "Adak")
col_leg <- c("gold", "#feb24c", "sienna", "tomato1", "orangered2", "#bd0026")


#to graph with lines between samples
s.class(pca_neutral$li, fac=pop(my_data_neutral), 
        col=col, #color of points. will retain lines between points
        clabel=0, #remove population labels
        cellipse=0, #remove ellipses; to add back in, make >=1
        cpoint=1,
        grid=FALSE, #otherwise will have light gray grid markers
        pch=19, #change point shapes
        axesell=TRUE)
# add legend
legend(x=0,y=60, legend = names_leg, col = col_leg, border = FALSE, bty = "n", cex = 0.9, pt.cex=1.5, title = "Sampling Site",pch=19)
# add eigenvalues plot as inset
add.scatter.eig(pca_neutral$eig[1:50],posi="bottomright", 3,2,1,ratio=.3)
# percent variation explained
eig.perc <- 100*pca_neutral$eig/sum(pca_neutral$eig)
head(eig.perc)




# Remove Outlier Individual -----------------------------------------------
install.packages("factoextra")
library(factoextra)

pca_neutral_coords <- get_pca_ind(pca_neutral)
pca_neutral_coords$coord

my_data_neutral2 <- read.genepop("batch_8_final_filtered_aligned_genepop_east_coastal_neutral_noUP45.gen")
X.n2 <- scaleGen(my_data_neutral2, NA.method="mean")
pca_neutral2 <- dudi.pca(X.n2,cent=FALSE,scale=FALSE,scannf=FALSE,nf=2) #nf = number to retain

#to graph with lines between samples
s.class(pca_neutral2$li, fac=pop(my_data_neutral2), 
        col=col, #color of points. will retain lines between points
        clabel=0, #remove population labels
        cellipse=0, #remove ellipses; to add back in, make >=1
        cpoint=1,
        grid=FALSE, #otherwise will have light gray grid markers
        pch=19, #change point shapes
        axesell=TRUE)
# add legend
legend(x=-70,y=75, legend = names_leg, col = col_leg, border = FALSE, bty = "n", cex = 0.9, pt.cex=1.5, title = "Sampling Site",pch=19)
# add eigenvalues plot as inset
add.scatter.eig(pca_neutral2$eig[1:50],posi="bottomleft", 3,2,1,ratio=.3)
# percent variation explained
eig.perc <- 100*pca_neutral2$eig/sum(pca_neutral2$eig)
head(eig.perc)
