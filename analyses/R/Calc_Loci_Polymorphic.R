############### Calculate % of loci that are polymorphic ################
#
# MF 5/7/2018
#
########################################################################



###################### WEST ###################
mydat <- read.delim("batch_8_final_filtered_aligned_GenotypesMat_west.txt", sep=" ", header=FALSE)


dim(mydat)
ncols(mydat)

poly <- 0
total <- 0
for(i in seq(2,ncol(mydat),1)){
  total = total + 1
  genotypes <- unique(mydat[,i])
  if(length(genotypes) > 2){
    poly = poly + 1
  }
}

poly / total


################### EAST ###################
mydat <- read.delim("batch_8_final_filtered_aligned_GenotypesMat_east.txt", sep=" ", header=FALSE)


dim(mydat)

poly <- 0
total <- 0
for(i in seq(2,ncol(mydat),1)){
  total = total + 1
  genotypes <- unique(mydat[,i])
  if(length(genotypes) > 2){
    poly = poly + 1
  }
}


unique(mydat[,2])
mydat[,3]


poly / total
