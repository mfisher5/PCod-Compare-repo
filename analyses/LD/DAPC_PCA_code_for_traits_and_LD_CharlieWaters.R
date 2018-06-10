#### Code for conducting and plotting DAPC's using optimal loci identified by Random Forest for each of my phenotypic traits

# Code for PCA's using corrected genotypes are at the bottom of the script


# April 13, 2016
# Charlie Waters

install.packages("adegenet")
install.packages("plotrix")
install.packages("ggplot2")
install.packages("genetics")
install.packages("LDheatmap")


library(adegenet)
library(plotrix)
library(ggplot2)
library(genetics)
library(LDheatmap)

setwd("C:/Users/Waters/Dropbox (MERLAB)/GWAS_Random_Forest/Random_Forest/Post_RF_analyses_using_optimal_loci")


################################################### DAPC using known groups

# Let's first run a DAPC with all individuals and all loci
data_all_loci <-read.genepop("Final_dataset_465indivs_9108loci_imputed_genepop.gen")
names(data_all_loci)
data_all_loci$pop

pop_98_Founders <- rep("P1",59)
pop_02int <- rep("F1 Wild",60)
pop_02seg <- rep("F1 Hatch", 55)
pop_06int <- rep("F2 INT",57)
pop_06seg <- rep("F2 SEG",53)
pop_10int <- rep("F3 INT",69)
pop_10seg <- rep("F3 SEG",59)
pop_14int <- rep("F4 INT",25)
pop_14seg <- rep("F4 SEG",28)

pop_groups <- as.factor(c(rep("P1",59),rep("F1 Wild",60),rep("F1 Hatch", 55),rep("F2 INT",57),rep("F2 SEG",53),rep("F3 INT",69),rep("F3 SEG",59),rep("F4 INT",25),rep("F4 SEG",28)))
pop_labels <- c(pop_98_Founders,pop_02int,pop_02seg,pop_06int,pop_06seg,pop_10int,pop_10seg,pop_14int,pop_14seg)
pop_cols <- c("black","dodgerblue","tomato","deepskyblue","red","blue","red4","deepskyblue4","sienna2")

dapc_all <- dapc(data_all_loci,data_all_loci$pop,n.pca=465,n.da=8) ##Retain all, then identify optimal number by optim.a.score
test_a_score <- optim.a.score(dapc_all)
dapc_all <- dapc(data_all_loci,data_all_loci$pop,n.pca=63,n.da=8) ##63 PC's is the optimal number

#2D plot
scatter(dapc_all,scree.da=FALSE,cellipse=0,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=1)
legend(x=4.7,y=3.8,bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)


#1D plot
par(mar=c(5,5,2,2))
scatter(dapc_all,1,1,scree.da=FALSE,cellipse=0,leg=FALSE,label=c("P1 Founders","F1 Wild","F1 Hatchery","F2 INT","F2 SEG","F3 INT","F3 SEG","F4 INT","F4 SEG"),posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=2,pch=c(18,15,15,16,16,17,17,15,15),solid=0.5)
legend("topright",legend=c("P1 Founders","F1 Wild","F1 Hatchery","F2 INT","F2 SEG","F3 INT","F3 SEG","F4 INT", "F4 SEG"),pch=c(15),col=pop_cols,cex=2)

#names(dapc_all)
#write.csv(dapc_all$var.contr,file="DAPC Variable Contributions.csv")

summary(dapc_all)  ####Gives summary of the DAPC - Assign.per.pop gives proportion of successful reassignments to original pops based on DF's
dapc_all$var  ### Proportion of variance conserved by the principal components
dapc_all$prior #### Numeric vector giving prior group probabilities
dapc_all$assign ## Posterior group assignment
dapc_all$posterior
dapc_all$eig[1]/sum(dapc_all$eig)  ### Variance explained by first discriminant function
dapc_all$eig[2]/sum(dapc_all$eig)  ### Variance explained by second discriminant function

# Identify structural loci of the DAPC
test_snpzip<-snpzip(data_all_loci,dapc_all,loading.plot=TRUE,method="median") 
test_snpzip
test_snpzip$FS[[1]][4]
test_snpzip$FS[1]
write.csv(test_snpzip$FS[[1]][4],file="Structural_Loci_Median_all_loci.csv")


###########################################################################################################################################################################################################
###########################################################################################################################################################################################################
###########################################################################################################################################################################################################
# Now conduct a DAPC using all trait-associated loci (for all traits combined, n=226 loci)
data_trait_loci <-read.genepop("Genepop_all_trait_associated_loci.gen")

dapc_traits <- dapc(data_trait_loci,data_trait_loci$pop,n.pca=226,n.da=8) ##Retain all, then identify optimal number by optim.a.score
test_a_score <- optim.a.score(dapc_traits)
dapc_traits <- dapc(data_trait_loci,data_trait_loci$pop,n.pca=118,n.da=8) ##118 PC's is the optimal number

#2D plot
scatter(dapc_traits,scree.da=FALSE,cellipse=0,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=1)
legend(x=4.4,y=4.1,bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)

dapc_traits$var  ### Proportion of variance conserved by the principal components 
dapc_traits$eig[1]/sum(dapc_traits$eig)  ### Variance explained by first discriminant function
dapc_traits$eig[2]/sum(dapc_traits$eig)  ### Variance explained by second discriminant function


####### Now compare with 226 random loci to see if trait loci separate pops more (possible evidence of selection)
data_random_226loci <-read.genepop("Genepop_random226_loci.gen")

dapc_random226 <- dapc(data_random_226loci,data_random_226loci$pop,n.pca=118,n.da=8) ##118 retain same number of PC's for direct comparison

#2D plot
scatter(dapc_random226,scree.da=FALSE,cellipse=0,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=1)
legend(x=3.7,y=3,bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)

dapc_random226$var  ### Proportion of variance conserved by the principal components 
dapc_random226$eig[1]/sum(dapc_random226$eig)  ### Variance explained by first discriminant function
dapc_random226$eig[2]/sum(dapc_random226$eig)  ### Variance explained by second discriminant function


###########################################################################################################################################################################################################
###########################################################################################################################################################################################################
###########################################################################################################################################################################################################

#Return Time

# read genepop file with all individuals and just the optimal loci for a trait (here, let's do return time)
return_time_optimal <-read.genepop("RF_optimal_return_genepop.gen")  #uncorrected genotypes

#or exclude founders to make comparable to PCA of corrected genotypes
return_time_optimal <-read.genepop("RF_optimal_return_genepop_no_founders.gen")  #uncorrected genotypes

summary(return_time_optimal) #no  missing data because we're using imputed data
names(return_time_optimal)
return_time_optimal$pop


#conduct a PCA on uncorrected genotypes using predictor loci for return time
return_uncorrected_pca <-dudi.pca(return_time_optimal,nf=10) #keep 10 axes
summary(return_uncorrected_pca) #First axis explains 7.04% variation; second axis explains 5.79% (if Founders included); If no founders, then axes explain 6.90 and 6.06%

#With founders
s.class(return_uncorrected_pca$li,pop(return_time_optimal),col=c("black","dodgerblue","tomato","deepskyblue","red","blue","red4","deepskyblue4","sienna2"),
        xax=1,yax=2,axesell=FALSE,clabel=1,cpoint=1,grid=FALSE,pch=c(18,15,15,16,16,17,17,15,15),
        label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"))

legend(x=3,y=8,bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),
       pch=c(18,15,15,16,16,17,17,15,15),col=c("black","dodgerblue","tomato","deepskyblue","red","blue","red4","deepskyblue4","sienna2"),cex=1.3)


#without founders
s.class(return_uncorrected_pca$li,pop(return_time_optimal),col=c("dodgerblue","tomato","deepskyblue","red","blue","red4","deepskyblue4","sienna2"),
        xax=1,yax=2,axesell=FALSE,clabel=1,cpoint=1,grid=FALSE,pch=c(15,15,16,16,17,17,15,15),
        label=c("F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"))

legend(x=2.8,y=6,bty='n',legend=c("F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),
       pch=c(15,15,16,16,17,17,15,15),col=c("dodgerblue","tomato","deepskyblue","red","blue","red4","deepskyblue4","sienna2"),cex=1.3)




# read in phenotypic data that has return time for each individual
pheno_data <- read.csv("Phenotypes_all_indivs_outliers_removed_missing_covariates_replaced.csv",header=TRUE)
return_time <- pheno_data$roza_day_of_year
hist(return_time,breaks=50)

# Use the BIC clustering criteria to identify the optimal number of clusters in the data
groups_return <- find.clusters(return_time_optimal) # retained 15 groups

#Now conduct a DAPC using those clusters
dapc_return <- dapc(return_time_optimal,groups_return$grp,n.pca=26,n.da=8) ##Retain all, then identify optimal number by optim.a.score
test_a_score <- optim.a.score(dapc_return)
dapc_return <- dapc(return_time_optimal,return_time_optimal$pop,n.pca=14,n.da=8) ##14 PC's is the optimal number

# DAPC using pre-defined populations
dapc_return <- dapc(return_time_optimal,return_time_optimal$pop,n.pca=26,n.da=8) ##Retain all, then identify optimal number by optim.a.score
test_a_score <- optim.a.score(dapc_return)
dapc_return <- dapc(return_time_optimal,return_time_optimal$pop,n.pca=26,n.da=8) ##26 PC's is the optimal number

#2D plot
scatter(dapc_return,scree.da=FALSE,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=1)
legend(x=2.5,y=3.5,bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)

#1D plot
scatter(dapc_return,1,1,scree.da=FALSE,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=0.5)
legend("topright",bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)
title("Return Time to Roza Dam")


dapc_return$var  ### Proportion of variance conserved by the principal components
dapc_return$eig[1]/sum(dapc_return$eig)  ### Variance explained by first discriminant function
dapc_return$eig[2]/sum(dapc_return$eig)  ### Variance explained by second discriminant function

# For return time, 100% variance is retained by PC's
# First DF explains 33.7% of retained variation
# Second DF explains 24.0% of retained variation



# Now try DAPC using 26 random loci

# read genepop file with all individuals and 26 random loci for return time
return_time_random <-read.genepop("Genepop_26random_loci.gen")  #uncorrected genotypes

# DAPC using pre-defined populations

dapc_return_random <- dapc(return_time_random,return_time_random$pop,n.pca=26,n.da=8) ##26 PC's is the optimal number

#2D plot
scatter(dapc_return_random,scree.da=FALSE,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=1)
legend(x=4.5,y=3.6,bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)

dapc_return_random$var  ### Proportion of variance conserved by the principal components
dapc_return_random$eig[1]/sum(dapc_return_random$eig)  ### Variance explained by first discriminant function
dapc_return_random$eig[2]/sum(dapc_return_random$eig)  ### Variance explained by second discriminant function






# The DAPC didn't give much separation, so now try a PCA but color the points on a scale from early to late arrival

#PCA on raw genotypes
pca_return <- dudi.pca(return_time_optimal,center=TRUE,scale=TRUE,nf=8)
return_coordinates <- pca_return$li  #extract the coordinates myself for plotting, because it gives me more control over plotting features



#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))
return_coordinates$Col <- rbPal(5)[as.numeric(cut(return_time,breaks = 5))] #adjust numbers to get different scales of colors
return_coordinates$day <- return_time
plot(return_coordinates$Axis1,return_coordinates$Axis2,pch=19,col=return_coordinates$Col)


#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))
pca_return2_coordinates$Col <- rbPal(10)[as.numeric(cut(return_optimal_corrected$roza_day_of_year,breaks = 10))]
pca_return2_coordinates$roza_day <- return_optimal_corrected$roza_day_of_year
plot(pca_return2_coordinates$PC1,pca_return2_coordinates$PC2,pch=19,col=pca_return2_coordinates$Col)

# Try each pop individually
return_pops <- seppop(return_time_optimal,pop=pop_groups)
return_1998 <- return_pops$P1
pca_return1998 <- dudi.pca(return_1998,center=TRUE,scale=TRUE)
summary(pca_return1998)
pca_return1998_coordinates <- pca_return1998$li  #extract the coordinates myself for plotting, because it gives me more control over plotting features

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))
pca_return1998_coordinates$Col <- rbPal(5)[as.numeric(cut(return_time[1:59],breaks = 5))] #adjust numbers to get different scales of colors
pca_return1998_coordinates$day <- return_time[1:59]
plot(pca_return1998_coordinates$Axis1,pca_return1998_coordinates$Axis2,pch=19,col=pca_return1998_coordinates$Col)


# Correspondence Analysis - This isn't a good technique because you lose a lot of within-pop info
return_genpop <- genind2genpop(return_time_optimal)
return_ca1 <- dudi.coa(tab(return_genpop),scannf=FALSE,nf=8)

barplot(return_ca1$eig,main="Correspondance Analysis Eigenvalues",
        col=heat.colors(length(return_ca1$eig)))

s.label(return_ca1$li, sub="CA 1-2",csub=2)
add.scatter.eig(return_ca1$eig,nf=8,xax=1,yax=2,posi="bottomright")






###########################################################################################################################################################################################################
###########################################################################################################################################################################################################
###########################################################################################################################################################################################################

#Spawn time

# read genepop file with all individuals and just the optimal loci for a trait (here, let's do spawn time)
spawn_time_optimal <-read.genepop("RF_optimal_spawn_genepop.gen")
summary(spawn_time_optimal)
names(spawn_time_optimal)
spawn_time_optimal$pop

# read in phenotypic data that has return time for each individual
pheno_data <- read.csv("Phenotypes_all_indivs_outliers_removed_missing_covariates_replaced.csv",header=TRUE)
spawn_time <- pheno_data$spawn_day_of_year
hist(spawn_time,breaks=50)

# DAPC using pre-defined populations
dapc_spawn <- dapc(spawn_time_optimal,spawn_time_optimal$pop,n.pca=68,n.da=8) ##Retain all, then identify optimal number by optim.a.score
test_a_score <- optim.a.score(dapc_spawn)
dapc_spawn <- dapc(spawn_time_optimal,spawn_time_optimal$pop,n.pca=31,n.da=8) ##31 PC's is the optimal number

#2D plot
scatter(dapc_spawn,scree.da=FALSE,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=1)
legend(x=3.4,y=3.6,bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)

#1D plot
scatter(dapc_spawn,1,1,scree.da=FALSE,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=0.5)
legend("topright",bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)
title("Spawn Timing")



dapc_spawn$var  ### Proportion of variance conserved by the principal components
dapc_spawn$eig[1]/sum(dapc_spawn$eig)  ### Variance explained by first discriminant function
dapc_spawn$eig[2]/sum(dapc_spawn$eig)  ### Variance explained by second discriminant function

# For spawn time, 78% % variance is retained by PC's
# First DF explains 44.2% of retained variation
# Second DF explains 16.8% of retained variation



# Now try DAPC using 68 random loci

# read genepop file with all individuals and 68 random loci for spawn time
spawn_time_random <-read.genepop("Genepop_68random_loci.gen")  #uncorrected genotypes

# DAPC using pre-defined populations

dapc_spawn_random <- dapc(spawn_time_random,spawn_time_random$pop,n.pca=31,n.da=8) ##31 PC's is the optimal number

#2D plot
scatter(dapc_spawn_random,scree.da=FALSE,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=1)
legend(x=2.7,y=3.2,bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)

dapc_spawn_random$var  ### Proportion of variance conserved by the principal components
dapc_spawn_random$eig[1]/sum(dapc_spawn_random$eig)  ### Variance explained by first discriminant function
dapc_spawn_random$eig[2]/sum(dapc_spawn_random$eig)  ### Variance explained by second discriminant function









# Use the BIC clustering criteria to identify the optimal number of clusters in the data
groups_spawn <- find.clusters(spawn_time_optimal) # keep all 68 PCs; identified 6 groups

#Now conduct a DAPC using those clusters
dapc_spawn <- dapc(spawn_time_optimal,groups_spawn$grp,n.pca=68,n.da=5) ##Retain all, then identify optimal number by optim.a.score
test_a_score <- optim.a.score(dapc_spawn)
dapc_spawn <- dapc(spawn_time_optimal,spawn_time_optimal$pop,n.pca=15,n.da=5) ##15 PC's is the optimal number

#2D plot
scatter(dapc_spawn)


# The DAPC didn't give much separation, so now try a PCA but color the points on a scale from early to late arrival

#PCA on raw genotypes
pca_spawn <- dudi.pca(spawn_time_optimal,center=TRUE,scale=TRUE,nf=8)
summary(pca_spawn)
spawn_coordinates <- pca_spawn$li  #extract the coordinates myself for plotting, because it gives me more control over plotting features
spawn_coordinates$spawn_day <- spawn_time

#omit individuals with missing spawn time 
spawn_coordinates_no_missing <- spawn_coordinates[which(spawn_coordinates$spawn_day!="NA"),]

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))
spawn_coordinates_no_missing$Col <- rbPal(5)[as.numeric(cut(spawn_coordinates_no_missing$spawn_day,breaks = 5))] #adjust numbers to get different scales of colors

plot(spawn_coordinates_no_missing$Axis2,spawn_coordinates_no_missing$Axis3,pch=19,col=spawn_coordinates_no_missing$Col)


##########################################################################################################################
# try dapc of spawn time when "pops" are the quartiles of spawn timing
# read genepop file with all individuals and just the optimal loci for a trait (here, let's do spawn time)
spawn_time_optimal <-read.genepop("RF_optimal_spawn_genepop_quartiles.gen")
names(spawn_time_optimal)
spawn_time_optimal$pop

first_quartile <- rep("Q1",131)
second_quartile <- rep("Q2",129)
third_quartile <- rep("Q3",87)
fourth_quartile <- rep("Q4",102)

quartile_labels <- c(first_quartile,second_quartile,third_quartile,fourth_quartile)
quartile_cols <- c("black","dodgerblue","tomato","purple")

#Now conduct a DAPC using those clusters
dapc_spawn <- dapc(spawn_time_optimal,spawn_time_optimal$pop,n.pca=68,n.da=3) ##Retain all, then identify optimal number by optim.a.score
test_a_score <- optim.a.score(dapc_spawn)
dapc_spawn <- dapc(spawn_time_optimal,spawn_time_optimal$pop,n.pca=44,n.da=3) ##44 PC's is the optimal number

#2D plot
scatter(dapc_spawn,col=quartile_cols,label=c("Q1","Q2","Q3","Q4"))

scatter(dapc_spawn,1,1,col=quartile_cols,label=c("Q1","Q2","Q3","Q4"))



###########################################################################################################################################################################################################
###########################################################################################################################################################################################################
###########################################################################################################################################################################################################

#Weight at Roza

# read genepop file with all individuals and just the optimal loci for a trait (here, let's do weight at roza)
weight_optimal <-read.genepop("RF_optimal_weight_genepop.gen")
names(weight_optimal)
weight_optimal$pop

# DAPC using pre-defined populations
dapc_weight <- dapc(weight_optimal,weight_optimal$pop,n.pca=37,n.da=8) ##Retain all, then identify optimal number by optim.a.score
test_a_score <- optim.a.score(dapc_weight)
dapc_weight <- dapc(weight_optimal,weight_optimal$pop,n.pca=37,n.da=8) ##37 PC's is the optimal number

#2D plot
scatter(dapc_weight,scree.da=FALSE,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=1)
legend(x=3.9,y=2.8,bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)

#1D plot
scatter(dapc_weight,1,1,scree.da=FALSE,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=1)
legend("topright",bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)
title("Weight at Roza Dam")



dapc_weight$var  ### Proportion of variance conserved by the principal components
dapc_weight$eig[1]/sum(dapc_weight$eig)  ### Variance explained by first discriminant function
dapc_weight$eig[2]/sum(dapc_weight$eig)  ### Variance explained by second discriminant function

# For weight, 100% variance is retained by PC's
# First DF explains 28.9% of retained variation
# Second DF explains 22.1% of retained variation


# read genepop file with all individuals and 37 random loci for weight
weight_time_random <-read.genepop("Genepop_37random_loci.gen")  #uncorrected genotypes

# DAPC using pre-defined populations

dapc_weight_random <- dapc(weight_time_random,weight_time_random$pop,n.pca=37,n.da=8) ##26 PC's is the optimal number

#2D plot
scatter(dapc_weight_random,scree.da=FALSE,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=1)
legend(x=2.2,y=3,bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)

dapc_weight_random$var  ### Proportion of variance conserved by the principal components
dapc_weight_random$eig[1]/sum(dapc_weight_random$eig)  ### Variance explained by first discriminant function
dapc_weight_random$eig[2]/sum(dapc_weight_random$eig)  ### Variance explained by second discriminant function







## Now do DAPC with "optimal loci" identified by a preliminary RF run with 
## no corrected genotypes but including covariates directly into RF

# read genepop file with all individuals and just the optimal loci for a trait (here, let's do weight at roza)
weight_optimal <-read.genepop("RF_preliminary_optimal_weight_genepop.gen")

# DAPC using pre-defined populations
dapc_weight <- dapc(weight_optimal,weight_optimal$pop,n.pca=52,n.da=8) ##Retain all, then identify optimal number by optim.a.score
test_a_score <- optim.a.score(dapc_weight)
dapc_weight <- dapc(weight_optimal,weight_optimal$pop,n.pca=52,n.da=8) ##52 PC's is the optimal number

#2D plot
scatter(dapc_weight,scree.da=FALSE,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=1)
legend(x=4.2,y=-0.8,bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)

#1D plot
scatter(dapc_weight,1,1,scree.da=FALSE,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=1)
legend("topright",bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)
title("Weight at Roza Dam")



dapc_weight$var  ### Proportion of variance conserved by the principal components
dapc_weight$eig[1]/sum(dapc_weight$eig)  ### Variance explained by first discriminant function
dapc_weight$eig[2]/sum(dapc_weight$eig)  ### Variance explained by second discriminant function

# For weight, 100% variance is retained by PC's
# First DF explains 28.9% of retained variation
# Second DF explains 22.1% of retained variation






###########################################################################################################################################################################################################
###########################################################################################################################################################################################################
###########################################################################################################################################################################################################

#Forklength at Roza

# read genepop file with all individuals and just the optimal loci for a trait (here, let's do forklength at roza)
forklength_optimal <-read.genepop("RF_optimal_forklength_genepop.gen") #uncorrected genotypes


#uncorrected genotypes without founders
forklength_optimal <-read.genepop("RF_optimal_forklength_genepop_no_founders.gen") 
names(forklength_optimal)
forklength_optimal$pop

#conduct a PCA on uncorrected genotypes using predictor loci for forklength
forklength_uncorrected_pca <-dudi.pca(forklength_optimal,nf=10) #keep 10 axes
summary(forklength_uncorrected_pca) #If no founders, then axes explain 4.78 and 3.84%


#without founders
s.class(forklength_uncorrected_pca$li,pop(forklength_optimal),col=c("dodgerblue","tomato","deepskyblue","red","blue","red4","deepskyblue4","sienna2"),
        xax=1,yax=2,axesell=FALSE,clabel=1,cpoint=1,grid=FALSE,pch=c(15,15,16,16,17,17,15,15),
        label=c("F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"))

legend(x=5,y=10,bty='n',legend=c("F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),
       pch=c(15,15,16,16,17,17,15,15),col=c("dodgerblue","tomato","deepskyblue","red","blue","red4","deepskyblue4","sienna2"),cex=1.3)




# DAPC using pre-defined populations
dapc_forklength <- dapc(forklength_optimal,forklength_optimal$pop,n.pca=44,n.da=8) ##Retain all, then identify optimal number by optim.a.score
test_a_score <- optim.a.score(dapc_forklength)
dapc_forklength <- dapc(forklength_optimal,forklength_optimal$pop,n.pca=25,n.da=8) ##25 PC's is the optimal number

#2D plot
scatter(dapc_forklength,scree.da=FALSE,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=1)
legend(x=2.5,y=3.3,bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)

#1D plot
scatter(dapc_forklength,1,1,scree.da=FALSE,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=0.5)
legend("topright",bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)
title("Forklength at Roza Dam")



dapc_forklength$var  ### Proportion of variance conserved by the principal components
dapc_forklength$eig[1]/sum(dapc_forklength$eig)  ### Variance explained by first discriminant function
dapc_forklength$eig[2]/sum(dapc_forklength$eig)  ### Variance explained by second discriminant function

# For forklength, 87.6% variance is retained by PC's
# First DF explains 40.8% of retained variation
# Second DF explains 20.3% of retained variation


# read genepop file with all individuals and 44 random loci for forklength
forklength_time_random <-read.genepop("Genepop_44random_loci.gen")  #uncorrected genotypes

# DAPC using pre-defined populations

dapc_forklength_random <- dapc(forklength_time_random,forklength_time_random$pop,n.pca=25,n.da=8) ##25 PC's is the optimal number

#2D plot
scatter(dapc_forklength_random,scree.da=FALSE,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=1)
legend(x=3.5,y=3,bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)

dapc_forklength_random$var  ### Proportion of variance conserved by the principal components
dapc_forklength_random$eig[1]/sum(dapc_forklength_random$eig)  ### Variance explained by first discriminant function
dapc_forklength_random$eig[2]/sum(dapc_forklength_random$eig)  ### Variance explained by second discriminant function





###########################################################################################################################################################################################################
###########################################################################################################################################################################################################
###########################################################################################################################################################################################################

#Daily growth coefficient

# read genepop file with all individuals and just the optimal loci for a trait (here, let's do daily growth coefficient)
DGC_optimal <-read.genepop("RF_optimal_DGC_genepop.gen")
names(DGC_optimal)
DGC_optimal$pop

# DAPC using pre-defined populations
dapc_DGC <- dapc(DGC_optimal,DGC_optimal$pop,n.pca=35,n.da=8) ##Retain all, then identify optimal number by optim.a.score
test_a_score <- optim.a.score(dapc_DGC)
dapc_DGC <- dapc(DGC_optimal,DGC_optimal$pop,n.pca=35,n.da=8) ##35 PC's is the optimal number

#2D plot
scatter(dapc_DGC,scree.da=FALSE,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=1)
legend(x=-5,y=3.4,bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)

#1D plot
scatter(dapc_DGC,1,1,scree.da=FALSE,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=0.5)
legend("topright",bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)
title("Daily Growth Coefficient")


dapc_DGC$var  ### Proportion of variance conserved by the principal components
dapc_DGC$eig[1]/sum(dapc_DGC$eig)  ### Variance explained by first discriminant function
dapc_DGC$eig[2]/sum(dapc_DGC$eig)  ### Variance explained by second discriminant function

# For DGC, 100% variance is retained by PC's
# First DF explains 25.1% of retained variation
# Second DF explains 17.5% of retained variation


# read genepop file with all individuals and 35 random loci for DGC
DGC_time_random <-read.genepop("Genepop_35random_loci.gen")  #uncorrected genotypes

# DAPC using pre-defined populations

dapc_DGC_random <- dapc(DGC_time_random,DGC_time_random$pop,n.pca=35,n.da=8) ##35 PC's is the optimal number

#2D plot
scatter(dapc_DGC_random,scree.da=FALSE,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=1)
legend(x=4.5,y=-0.5,bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)

dapc_DGC_random$var  ### Proportion of variance conserved by the principal components
dapc_DGC_random$eig[1]/sum(dapc_DGC_random$eig)  ### Variance explained by first discriminant function
dapc_DGC_random$eig[2]/sum(dapc_DGC_random$eig)  ### Variance explained by second discriminant function




###########################################################################################################################################################################################################
###########################################################################################################################################################################################################
###########################################################################################################################################################################################################

#Age at maturity (using optimal loci identified by RF when genotypes but not phenotypes were corrected and having a balanced sampling design in RF)

# read genepop file with all individuals and just the optimal loci for a trait (here, let's do age at maturity)
age_optimal <-read.genepop("RF_optimal_age_genepop.gen")
names(age_optimal)
age_optimal$pop

# DAPC using pre-defined populations
dapc_age <- dapc(age_optimal,age_optimal$pop,n.pca=30,n.da=8) ##Retain all, then identify optimal number by optim.a.score
test_a_score <- optim.a.score(dapc_age)
dapc_age <- dapc(age_optimal,age_optimal$pop,n.pca=24,n.da=8) ##24 PC's is the optimal number

#2D plot
scatter(dapc_age,scree.da=FALSE,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=1)
legend(x=-5.5,y=3,bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)

#1D plot
scatter(dapc_age,1,1,scree.da=FALSE,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=0.5)
legend("topright",bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)
title("Age at Maturity")


dapc_age$var  ### Proportion of variance conserved by the principal components
dapc_age$eig[1]/sum(dapc_age$eig)  ### Variance explained by first discriminant function
dapc_age$eig[2]/sum(dapc_age$eig)  ### Variance explained by second discriminant function

# For age, 95.0% variance is retained by PC's
# First DF explains 28.7% of retained variation
# Second DF explains 21.3% of retained variation


# read genepop file with all individuals and 30 random loci for age at maturity
age_time_random <-read.genepop("Genepop_30random_loci.gen")  #uncorrected genotypes

# DAPC using pre-defined populations

dapc_age_random <- dapc(age_time_random,age_time_random$pop,n.pca=24,n.da=8) ##24 PC's is the optimal number

#2D plot
scatter(dapc_age_random,scree.da=FALSE,leg=FALSE,label=c("P1","F1W","F1H","F2I","F2S","F3I","F3S","F4I","F4S"),
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,clabel=1,pch=c(18,15,15,16,16,17,17,15,15),solid=1)
legend(x=-5.1,y=2.7,bty='n',legend=c("P1: Founders","F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),pch=c(18,15,15,16,16,17,17,15,15),col=pop_cols,cex=1.3)

dapc_age_random$var  ### Proportion of variance conserved by the principal components
dapc_age_random$eig[1]/sum(dapc_age_random$eig)  ### Variance explained by first discriminant function
dapc_age_random$eig[2]/sum(dapc_age_random$eig)  ### Variance explained by second discriminant function


###############################################################################
# Now try PCA on corrected genotypes and phenotypes (since RF was conducted on those) and including only optimal loci using PCA

# Return Time
return_optimal_corrected <- read.csv("return_time_optimal_loci_corrected_phenos_genos.csv",header=TRUE,row.names=1)
return_genos <- return_optimal_corrected[,3:28]
pca_return2 <- prcomp(return_genos)
summary(pca_return2)  #PC 1 explains 9.9% of variation; PC 2 is 9.1%

pca_return2_coordinates <- as.data.frame(pca_return2$x)

return_time_pops <- as.factor(c(rep("F1W",54),rep("F1H",50),rep("F2I",53),
                                rep("F2S",49),rep("F3I",65),rep("F3S",59),
                                rep("F4I",25),rep("F4S",28)))

return_time_cols <- c(rep("dodgerblue",54),rep("tomato",50),rep("deepskyblue",53),
                                rep("red",49),rep("blue",65),rep("red4",59),
                                rep("deepskyblue4",25),rep("sienna2",28))

#For s.class, the pop colors must match the order of the populations when they are listed as factors (so the order of the factor levels)
s.class(pca_return2_coordinates,return_time_pops,col=c("tomato","dodgerblue","deepskyblue","red","blue","red4","deepskyblue4","sienna2"),
        xax=1,yax=2,axesell=FALSE,clabel=1,cpoint=1,grid=FALSE,pch=c(15,15,16,16,17,17,15,15))

legend(x=2.25,y=2.5,bty='n',legend=c("F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),
       pch=c(15,15,16,16,17,17,15,15),col=c("dodgerblue","tomato","deepskyblue","red","blue","red4","deepskyblue4","sienna2"),cex=1.3)


#Check to make sure s.class plot is accurate
plot(pca_return2_coordinates$PC1,pca_return2_coordinates$PC2,col=return_time_cols,pch=19)



##################################################################################################################################
# Forklength at Roza
length_optimal_corrected <- read.csv("Forklength_optimal_loci_corrected_all_covariates.csv",header=TRUE,row.names=1)
length_genos <- length_optimal_corrected[,3:46]
pca_length <- prcomp(length_genos)
summary(pca_length)  # PC1 explains 7.5% of the variation; PC2 explains 6.9%
pca_length_coordinates <- as.data.frame(pca_length$x)

forklength_pops <- as.factor(c(rep("F1W",52),rep("F1H",50),rep("F2I",53),
                                rep("F2S",49),rep("F3I",65),rep("F3S",58),
                                rep("F4I",25),rep("F4S",28)))

forklength_cols <- c(rep("dodgerblue",52),rep("tomato",50),rep("deepskyblue",53),
                      rep("red",49),rep("blue",65),rep("red4",58),
                      rep("deepskyblue4",25),rep("sienna2",28))

#For s.class, the pop colors must match the order of the populations when they are listed as factors (so the order of the factor levels)
s.class(pca_length_coordinates,forklength_pops,col=c("tomato","dodgerblue","deepskyblue","red","blue","red4","deepskyblue4","sienna2"),
        xax=1,yax=2,axesell=FALSE,clabel=1,cpoint=1,grid=FALSE,pch=c(15,15,16,16,17,17,15,15))

legend(x=2.05,y=2.25,bty='n',legend=c("F1W: F1 Wild","F1H: F1 Hatchery","F2I: F2 INT","F2S: F2 SEG","F3I: F3 INT","F3S: F3 SEG","F4I: F4 INT","F4S: F4 SEG"),
       pch=c(15,15,16,16,17,17,15,15),col=c("dodgerblue","tomato","deepskyblue","red","blue","red4","deepskyblue4","sienna2"),cex=1.3)


#Check to make sure s.class plot is accurate
plot(pca_length_coordinates$PC1,pca_length_coordinates$PC2,col=forklength_cols,pch=19)

















# Calculate linkage disequilibrium for predictor loci

#Return Time

# read file with all individuals and just the optimal loci for a trait (here, let's do return time)
# genepop format was replaced with AA, AT, and TT to conform with the genetics package format (i.e. 0101 is AA, 0102 is AT, and 0202 is TT)

#loci for return time are ordered by map position (unmapped loci are at the end)
return_time_optimal <- read.table("optimal_loci_return_time_genetics_format.txt",stringsAsFactors=FALSE) #need the string option so that each column can be a vector

return_time_genos <- return_time_optimal[2:466,2:27]

#create dummy data frame
return_time_geno_all <- data.frame(matrix(NA, nrow = 465, ncol = 26)) #make the same dimensions as the "genos" file in the previous line

for (i in 1:ncol(return_time_genos)) {
  return_time_locus <- unname(unlist(return_time_genos[,i]))
  return_time_locus_geno <- genotype(return_time_locus)
  return_time_geno_all[,i] <- return_time_locus_geno
  
}

return_time_geno_all  

# Rows for different pops are as follows:
# 1998 founders: 
# 2002 INT:
# 2002 SEG:
# 2006 INT:
# 2006 SEG:
# 2010 INT: 285-353
# 2010 SEG: 354-412
# 2014 INT:
# 2014 SEG:

colnames(return_time_geno_all) <- return_time_optimal[1,2:27]
return_time_heatmap <- LDheatmap(return_time_geno_all,SNP.name=unname(unlist(return_time_optimal[1,2:27])),
                                 add.map=FALSE,color=heat.colors(20),title="Pairwise LD-Return Time Predictor Loci")








