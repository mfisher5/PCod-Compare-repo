############# Sliding Window Analysis #########
#
# Adopted from Charlie Waters' script "Code_Kernel_smoothing_4_generations"
# Original function written by Kotaro Ono
#
# MF edited for Pacific cod 3/14/2018
#
###############################################



# Install Packages --------------------------------------------------------

install.packages("TTR")
library(TTR) 
library(zoo) # don't need to install; will install with TTR

install.packages("dplyr") #always install this last!
library(dplyr)


# Load Data ---------------------------------------------------------------

## set working directory
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/SlidingWindow")

## read in files with loci positions and fst
east <- read.delim("batch_8_SWA_input_east_sorted.txt",header=TRUE)
head(east)
dim(east)

west <- read.delim("West/batch_8_SWA_input_west_sorted.txt",sep="\t", header=TRUE)
head(west)
dim(west)




# Match Data Sets by Locus ------------------------------------------------

## join data files into a single frame. inner join to match loci that exist in both groups
align_data <- inner_join(east,west,by="Locus")
dim(align_data)
head(align_data)

## split data files back into east and west
east_matched <- select(align_data, c("Locus", "fst.x", "chromosome.x", "position.x"))
head(east_matched)
colnames(east_matched) <- c("Locus", "fst", "chromosome", "position")

west_matched <- select(align_data, c("Locus", "fst.y", "chromosome.y", "position.y"))
head(west_matched)
colnames(west_matched) <- c("Locus", "fst", "chromosome", "position")
dim(east_matched)
dim(west_matched)



# Write out Data ----------------------------------------------------------
## write out to file to save dataframes
write.table(east_matched, "batch_8_SWA_input_east_matched.txt", quote=FALSE, sep="\t")
write.table(west_matched, "batch_8_SWA_input_west_matched.txt", quote=FALSE, sep="\t")


# Create SLA Function (no plotting) -----------------------------------------------------

## note how there are defaults suppiled to each argument. these defaults are based on Charlie's script.

plotting_reg_interval = function(Dat=mydata, Sigma_sliding_window=3, Nb_bootstrap=1000000, which.chromosome.analysis="all", which.chromosome.plot="All",division=100,name_output="SLA_output", path_output=".") {
  if(which.chromosome.analysis!="all") Data.analysis = subset(Dat, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis=="all") Data.analysis = Dat
  
  N_data = nrow(Data.analysis)
  Nb_chromosome = unique(Data.analysis$chromosome) 
  
  ## initialize parameters
  Expected_MA = c()
  Bootstrap_mean = c()
  Bootstrap_sd = c()
  CI_all_chromosomes=c()
  Position_MA = c()
  to_export=c()
  to_export=as.data.frame(to_export)
  
  ## for each chromosome:
  for (ijk in Nb_chromosome)
  {
    Data.part = subset(Dat, Dat$chromosome %in% ijk)
    Size_window = 3*Sigma_sliding_window  # sliding window is 3 times the value of sigma in each direction
    Beg_position_cM=min(Data.part$position)+Size_window ## calculations cannot start at 0, because it would be biased no values between -6 and 0cM for example
    End_position_cM=max(Data.part$position)-Size_window
    Beg_position = which(Data.part$position >=Beg_position_cM)[1] ### start at the first position as calculated before or at the first position where there is a data point (as sometimes, a large portion of the chromosome is made of duplicated loci --> not kept for population studies)
    End_position = which(Data.part$position >= (max(End_position_cM)))[1]
    Positions_LG=seq(from=Beg_position_cM,to=End_position_cM,length.out=division)  ## these position will be equidistant
    MA_calc = function(x, do.boot=FALSE) ## MA=Moving Average
    {
      if(do.boot==TRUE) Fst_data = Dat$fst[sample(1:nrow(Dat), size= nrow(Data.part), replace=TRUE)]
      if(do.boot==FALSE) Fst_data = Data.part$fst
      MA_windows = which(Data.part$position >= (Positions_LG[x]-Size_window) & Data.part$position <= (Positions_LG[x]+Size_window))
      MA_weight = exp(-(Data.part$position[MA_windows]-Positions_LG[x])^2/(2*Sigma_sliding_window^2))
      MA_val = sum(Fst_data[MA_windows]*MA_weight)/sum(MA_weight)
    }  
    
    # calculate the MAs_exp for all positions
    MAs_exp = sapply(1:division, function(x) MA_calc(x, do.boot=FALSE))[1:division]      
    
    # calculate the MAs_bootstrap for all positions  
    MAs_bootstrap = sapply(1:Nb_bootstrap, function(x) sapply(1:division, function(y) MA_calc(y, do.boot=TRUE))[1:division])
    #Calculate Confidence Interval
    CI_all=t(apply(MAs_bootstrap,1,function(x) quantile(x,c(0.5,0.025,0.975),na.rm=TRUE)))
    
    
    # Mean_boot = apply(MAs_bootstrap, 1, mean)
    Sd_boot = apply(MAs_bootstrap, 1, sd)
    
    Expected_MA = c(Expected_MA, MAs_exp)
    CI_all_chromosomes=rbind(CI_all_chromosomes,CI_all)
    Bootstrap_sd=c(Bootstrap_sd,Sd_boot)
    
    min_lower95=min(CI_all,na.rm=TRUE)
    max_upper95=max(CI_all,na.rm=TRUE)
    min_MA=min(MAs_exp,na.rm=TRUE)
    max_MA=max(MAs_exp,na.rm=TRUE)
    min_lim=min(c(min_lower95,max_upper95,min_MA,max_MA),na.rm=TRUE)-0.02
    max_lim=max(c(min_lower95,max_upper95,min_MA,max_MA),na.rm=TRUE)+0.02
    where.NA = which(is.na(CI_all[,2]))
    #plot(Positions_LG, MAs_exp, type="l", ylim=c(min_lim,max_lim), lwd=3, col="red", main = paste(ijk))
    #lines(Positions_LG, CI_all[,1], col=1, lwd=2)
    #lines(Positions_LG, CI_all[,2], col=rgb(0.3,0.1,0.4, 0.5))
    #lines(Positions_LG, CI_all[,3], col=rgb(0.3,0.1,0.4, 0.5))
    #if(length(where.NA)!=0) polygon(c(Positions_LG[-where.NA], rev(Positions_LG[-where.NA])), c(CI_all[,2][-where.NA], rev(CI_all[,3][-where.NA])), col=rgb(0.3,0.1,0.4, 0.5))
    #if(length(where.NA)==0) polygon(c(Positions_LG, rev(Positions_LG)), c(CI_all[,2], rev(CI_all[,3])), col=rgb(0.3,0.1,0.4, 0.5))
    positions_to_add=Positions_LG
    chromosome=matrix(ijk,length(positions_to_add),1)
    to_export_temp=cbind(chromosome,positions_to_add,MAs_exp,CI_all)
    to_export_temp=as.data.frame(to_export_temp)
    to_export=rbind(to_export,to_export_temp)
    print("Done with Chromosome:")
    print(ijk)
  }  
  
  
  (P_value = pnorm(Expected_MA, CI_all_chromosomes[,1], Bootstrap_sd))
  print (P_value)
  to_export=cbind(to_export,P_value)
  output<-paste(path_output,'/',name_output,'_kernel_smoothing_',Nb_bootstrap,'_bootstraps_sigma_',Sigma_sliding_window,'_div',division,'.txt', sep='')
  colnames(to_export)=c("chromosome","position","Fst/Fct","Mean_boostrap","lower_95","upper_95","pvalue")
  print(head(to_export))
  write.table(to_export,file=output,quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
}



# Run Sliding Window Analysis ---------------------------------------------
plotting_reg_interval(Dat = west_matched, Sigma_sliding_window=250000, Nb_bootstrap=100000, which.chromosome.analysis="all",which.chromosome.plot="all",division=150,name_output="West/batch_8_final_filtered_west_matched",path_output=".")
 




