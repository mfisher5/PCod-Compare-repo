setwd("C:/Users/Waters/Dropbox (MERLAB)/CESRF_2014_adult_sequencing/final_analysis/sliding_window_analysis")

install.packages("TTR")
install.packages("zoo")

library(TTR)      # to use the moving average function to speed-up the work
library(zoo)

###############################################################################################################################################################################################

# 2002INT sigma 3 (window=18cM); 1,000,000 reps

par(mfrow=c(6,6), mar=c(2,2,2,2)) ### change to correspond to your data
plotting_reg_interval = function(Sigma_sliding_window=3, Nb_bootstrap=1000000, which.chromosome.analysis="all", which.chromosome.plot="All",division=100,name_output="2002INT") ### these parameters should be specified. If they are not, the default will be used
{
  
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
    plot(Positions_LG, MAs_exp, type="l", ylim=c(min_lim,max_lim), lwd=3, col="red", main = paste(ijk))
    lines(Positions_LG, CI_all[,1], col=1, lwd=2)
    lines(Positions_LG, CI_all[,2], col=rgb(0.3,0.1,0.4, 0.5))
    lines(Positions_LG, CI_all[,3], col=rgb(0.3,0.1,0.4, 0.5))
    if(length(where.NA)!=0) polygon(c(Positions_LG[-where.NA], rev(Positions_LG[-where.NA])), c(CI_all[,2][-where.NA], rev(CI_all[,3][-where.NA])), col=rgb(0.3,0.1,0.4, 0.5))
    if(length(where.NA)==0) polygon(c(Positions_LG, rev(Positions_LG)), c(CI_all[,2], rev(CI_all[,3])), col=rgb(0.3,0.1,0.4, 0.5))
    positions_to_add=Positions_LG
    chromosome=matrix(ijk,length(positions_to_add),1)
    to_export_temp=cbind(chromosome,positions_to_add,MAs_exp,CI_all)
    to_export_temp=as.data.frame(to_export_temp)
    to_export=rbind(to_export,to_export_temp)
  }  
  
  
  (P_value = pnorm(Expected_MA, CI_all_chromosomes[,1], Bootstrap_sd))
  print (P_value)
  to_export=cbind(to_export,P_value)
  output<-paste('C:/Users/Waters/Dropbox (MERLAB)/CESRF_2014_adult_sequencing/final_analysis/sliding_window_analysis/',name_output,'_kernel_smoothing_',Nb_bootstrap,'_bootstraps_sigma_',Sigma_sliding_window,'_div',division,'.txt', sep='')
  colnames(to_export)=c("chromosome","position","Fst/Fct","Mean_boostrap","lower_95","upper_95","pvalue")
  print(head(to_export))
  write.table(to_export,file=output,quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
}



#windows()
Dat = read.table("sliding_window_input_2002int.txt", header=TRUE, sep="\t")
colnames(Dat) = c("Locus","fst", "chromosome", "position")
### Re-order the dataset so that for each chromosome, I have data ordered by the closest to furthest cm from the telomere
Dat = Dat[order(Dat$chromosome, Dat$position),]


par(mfrow=c(6,6), mar=c(2,2,2,2))
plotting_reg_interval(Sigma_sliding_window=3, Nb_bootstrap=1000000, which.chromosome.analysis="all",division=100,name_output="2002INT")


###############################################################################################################################################################################################

# 2002SEG sigma 3 (window=18cM); 1,000,000 reps

par(mfrow=c(6,6), mar=c(2,2,2,2)) ### change to correspond to your data
plotting_reg_interval = function(Sigma_sliding_window=3, Nb_bootstrap=1000000, which.chromosome.analysis="all", which.chromosome.plot="All",division=100,name_output="2002SEG") ### these parameters should be specified. If they are not, the default will be used
{
  
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
    plot(Positions_LG, MAs_exp, type="l", ylim=c(min_lim,max_lim), lwd=3, col="red", main = paste(ijk))
    lines(Positions_LG, CI_all[,1], col=1, lwd=2)
    lines(Positions_LG, CI_all[,2], col=rgb(0.3,0.1,0.4, 0.5))
    lines(Positions_LG, CI_all[,3], col=rgb(0.3,0.1,0.4, 0.5))
    if(length(where.NA)!=0) polygon(c(Positions_LG[-where.NA], rev(Positions_LG[-where.NA])), c(CI_all[,2][-where.NA], rev(CI_all[,3][-where.NA])), col=rgb(0.3,0.1,0.4, 0.5))
    if(length(where.NA)==0) polygon(c(Positions_LG, rev(Positions_LG)), c(CI_all[,2], rev(CI_all[,3])), col=rgb(0.3,0.1,0.4, 0.5))
    positions_to_add=Positions_LG
    chromosome=matrix(ijk,length(positions_to_add),1)
    to_export_temp=cbind(chromosome,positions_to_add,MAs_exp,CI_all)
    to_export_temp=as.data.frame(to_export_temp)
    to_export=rbind(to_export,to_export_temp)
  }  
  
  
  (P_value = pnorm(Expected_MA, CI_all_chromosomes[,1], Bootstrap_sd))
  print (P_value)
  to_export=cbind(to_export,P_value)
  output<-paste('C:/Users/Waters/Dropbox (MERLAB)/CESRF_2014_adult_sequencing/final_analysis/sliding_window_analysis/',name_output,'_kernel_smoothing_',Nb_bootstrap,'_bootstraps_sigma_',Sigma_sliding_window,'_div',division,'.txt', sep='')
  colnames(to_export)=c("chromosome","position","Fst/Fct","Mean_boostrap","lower_95","upper_95","pvalue")
  print(head(to_export))
  write.table(to_export,file=output,quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
}



#windows()
Dat = read.table("sliding_window_input_2002seg.txt", header=TRUE, sep="\t")
colnames(Dat) = c("Locus","fst", "chromosome", "position")
### Re-order the dataset so that for each chromosome, I have data ordered by the closest to furthest cm from the telomere
Dat = Dat[order(Dat$chromosome, Dat$position),]


par(mfrow=c(6,6), mar=c(2,2,2,2))
plotting_reg_interval(Sigma_sliding_window=3, Nb_bootstrap=1000000, which.chromosome.analysis="all",division=100,name_output="2002SEG")


###############################################################################################################################################################################################

# 2006INT sigma 3 (window=18cM); 1,000,000 reps

par(mfrow=c(6,6), mar=c(2,2,2,2)) ### change to correspond to your data
plotting_reg_interval = function(Sigma_sliding_window=3, Nb_bootstrap=1000000, which.chromosome.analysis="all", which.chromosome.plot="All",division=100,name_output="2006INT") ### these parameters should be specified. If they are not, the default will be used
{
  
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
    plot(Positions_LG, MAs_exp, type="l", ylim=c(min_lim,max_lim), lwd=3, col="red", main = paste(ijk))
    lines(Positions_LG, CI_all[,1], col=1, lwd=2)
    lines(Positions_LG, CI_all[,2], col=rgb(0.3,0.1,0.4, 0.5))
    lines(Positions_LG, CI_all[,3], col=rgb(0.3,0.1,0.4, 0.5))
    if(length(where.NA)!=0) polygon(c(Positions_LG[-where.NA], rev(Positions_LG[-where.NA])), c(CI_all[,2][-where.NA], rev(CI_all[,3][-where.NA])), col=rgb(0.3,0.1,0.4, 0.5))
    if(length(where.NA)==0) polygon(c(Positions_LG, rev(Positions_LG)), c(CI_all[,2], rev(CI_all[,3])), col=rgb(0.3,0.1,0.4, 0.5))
    positions_to_add=Positions_LG
    chromosome=matrix(ijk,length(positions_to_add),1)
    to_export_temp=cbind(chromosome,positions_to_add,MAs_exp,CI_all)
    to_export_temp=as.data.frame(to_export_temp)
    to_export=rbind(to_export,to_export_temp)
  }  
  
  
  (P_value = pnorm(Expected_MA, CI_all_chromosomes[,1], Bootstrap_sd))
  print (P_value)
  to_export=cbind(to_export,P_value)
  output<-paste('C:/Users/Waters/Dropbox (MERLAB)/CESRF_2014_adult_sequencing/final_analysis/sliding_window_analysis/',name_output,'_kernel_smoothing_',Nb_bootstrap,'_bootstraps_sigma_',Sigma_sliding_window,'_div',division,'.txt', sep='')
  colnames(to_export)=c("chromosome","position","Fst/Fct","Mean_boostrap","lower_95","upper_95","pvalue")
  print(head(to_export))
  write.table(to_export,file=output,quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
}



#windows()
Dat = read.table("sliding_window_input_2006int.txt", header=TRUE, sep="\t")
colnames(Dat) = c("Locus","fst", "chromosome", "position")
### Re-order the dataset so that for each chromosome, I have data ordered by the closest to furthest cm from the telomere
Dat = Dat[order(Dat$chromosome, Dat$position),]


par(mfrow=c(6,6), mar=c(2,2,2,2))
plotting_reg_interval(Sigma_sliding_window=3, Nb_bootstrap=1000000, which.chromosome.analysis="all",division=100,name_output="2006INT")


###############################################################################################################################################################################################

# 2006SEG sigma 3 (window=18cM); 1,000,000 reps

par(mfrow=c(6,6), mar=c(2,2,2,2)) ### change to correspond to your data
plotting_reg_interval = function(Sigma_sliding_window=3, Nb_bootstrap=1000000, which.chromosome.analysis="all", which.chromosome.plot="All",division=100,name_output="2006SEG") ### these parameters should be specified. If they are not, the default will be used
{
  
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
    plot(Positions_LG, MAs_exp, type="l", ylim=c(min_lim,max_lim), lwd=3, col="red", main = paste(ijk))
    lines(Positions_LG, CI_all[,1], col=1, lwd=2)
    lines(Positions_LG, CI_all[,2], col=rgb(0.3,0.1,0.4, 0.5))
    lines(Positions_LG, CI_all[,3], col=rgb(0.3,0.1,0.4, 0.5))
    if(length(where.NA)!=0) polygon(c(Positions_LG[-where.NA], rev(Positions_LG[-where.NA])), c(CI_all[,2][-where.NA], rev(CI_all[,3][-where.NA])), col=rgb(0.3,0.1,0.4, 0.5))
    if(length(where.NA)==0) polygon(c(Positions_LG, rev(Positions_LG)), c(CI_all[,2], rev(CI_all[,3])), col=rgb(0.3,0.1,0.4, 0.5))
    positions_to_add=Positions_LG
    chromosome=matrix(ijk,length(positions_to_add),1)
    to_export_temp=cbind(chromosome,positions_to_add,MAs_exp,CI_all)
    to_export_temp=as.data.frame(to_export_temp)
    to_export=rbind(to_export,to_export_temp)
  }  
  
  
  (P_value = pnorm(Expected_MA, CI_all_chromosomes[,1], Bootstrap_sd))
  print (P_value)
  to_export=cbind(to_export,P_value)
  output<-paste('C:/Users/Waters/Dropbox (MERLAB)/CESRF_2014_adult_sequencing/final_analysis/sliding_window_analysis/',name_output,'_kernel_smoothing_',Nb_bootstrap,'_bootstraps_sigma_',Sigma_sliding_window,'_div',division,'.txt', sep='')
  colnames(to_export)=c("chromosome","position","Fst/Fct","Mean_boostrap","lower_95","upper_95","pvalue")
  print(head(to_export))
  write.table(to_export,file=output,quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
}



#windows()
Dat = read.table("sliding_window_input_2006seg.txt", header=TRUE, sep="\t")
colnames(Dat) = c("Locus","fst", "chromosome", "position")
### Re-order the dataset so that for each chromosome, I have data ordered by the closest to furthest cm from the telomere
Dat = Dat[order(Dat$chromosome, Dat$position),]


par(mfrow=c(6,6), mar=c(2,2,2,2))
plotting_reg_interval(Sigma_sliding_window=3, Nb_bootstrap=1000000, which.chromosome.analysis="all",division=100,name_output="2006SEG")


###############################################################################################################################################################################################

# 2010INT sigma 3 (window=18cM); 1,000,000 reps

par(mfrow=c(6,6), mar=c(2,2,2,2)) ### change to correspond to your data
plotting_reg_interval = function(Sigma_sliding_window=3, Nb_bootstrap=1000000, which.chromosome.analysis="all", which.chromosome.plot="All",division=100,name_output="2010INT") ### these parameters should be specified. If they are not, the default will be used
{
  
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
    plot(Positions_LG, MAs_exp, type="l", ylim=c(min_lim,max_lim), lwd=3, col="red", main = paste(ijk))
    lines(Positions_LG, CI_all[,1], col=1, lwd=2)
    lines(Positions_LG, CI_all[,2], col=rgb(0.3,0.1,0.4, 0.5))
    lines(Positions_LG, CI_all[,3], col=rgb(0.3,0.1,0.4, 0.5))
    if(length(where.NA)!=0) polygon(c(Positions_LG[-where.NA], rev(Positions_LG[-where.NA])), c(CI_all[,2][-where.NA], rev(CI_all[,3][-where.NA])), col=rgb(0.3,0.1,0.4, 0.5))
    if(length(where.NA)==0) polygon(c(Positions_LG, rev(Positions_LG)), c(CI_all[,2], rev(CI_all[,3])), col=rgb(0.3,0.1,0.4, 0.5))
    positions_to_add=Positions_LG
    chromosome=matrix(ijk,length(positions_to_add),1)
    to_export_temp=cbind(chromosome,positions_to_add,MAs_exp,CI_all)
    to_export_temp=as.data.frame(to_export_temp)
    to_export=rbind(to_export,to_export_temp)
  }  
  
  
  (P_value = pnorm(Expected_MA, CI_all_chromosomes[,1], Bootstrap_sd))
  print (P_value)
  to_export=cbind(to_export,P_value)
  output<-paste('C:/Users/Waters/Dropbox (MERLAB)/CESRF_2014_adult_sequencing/final_analysis/sliding_window_analysis/',name_output,'_kernel_smoothing_',Nb_bootstrap,'_bootstraps_sigma_',Sigma_sliding_window,'_div',division,'.txt', sep='')
  colnames(to_export)=c("chromosome","position","Fst/Fct","Mean_boostrap","lower_95","upper_95","pvalue")
  print(head(to_export))
  write.table(to_export,file=output,quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
}



#windows()
Dat = read.table("sliding_window_input_2010int.txt", header=TRUE, sep="\t")
colnames(Dat) = c("Locus","fst", "chromosome", "position")
### Re-order the dataset so that for each chromosome, I have data ordered by the closest to furthest cm from the telomere
Dat = Dat[order(Dat$chromosome, Dat$position),]


par(mfrow=c(6,6), mar=c(2,2,2,2))
plotting_reg_interval(Sigma_sliding_window=3, Nb_bootstrap=1000000, which.chromosome.analysis="all",division=100,name_output="2010INT")


###############################################################################################################################################################################################

# 2010SEG sigma 3 (window=18cM); 1,000,000 reps

par(mfrow=c(6,6), mar=c(2,2,2,2)) ### change to correspond to your data
plotting_reg_interval = function(Sigma_sliding_window=3, Nb_bootstrap=1000000, which.chromosome.analysis="all", which.chromosome.plot="All",division=100,name_output="2010SEG") ### these parameters should be specified. If they are not, the default will be used
{
  
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
    plot(Positions_LG, MAs_exp, type="l", ylim=c(min_lim,max_lim), lwd=3, col="red", main = paste(ijk))
    lines(Positions_LG, CI_all[,1], col=1, lwd=2)
    lines(Positions_LG, CI_all[,2], col=rgb(0.3,0.1,0.4, 0.5))
    lines(Positions_LG, CI_all[,3], col=rgb(0.3,0.1,0.4, 0.5))
    if(length(where.NA)!=0) polygon(c(Positions_LG[-where.NA], rev(Positions_LG[-where.NA])), c(CI_all[,2][-where.NA], rev(CI_all[,3][-where.NA])), col=rgb(0.3,0.1,0.4, 0.5))
    if(length(where.NA)==0) polygon(c(Positions_LG, rev(Positions_LG)), c(CI_all[,2], rev(CI_all[,3])), col=rgb(0.3,0.1,0.4, 0.5))
    positions_to_add=Positions_LG
    chromosome=matrix(ijk,length(positions_to_add),1)
    to_export_temp=cbind(chromosome,positions_to_add,MAs_exp,CI_all)
    to_export_temp=as.data.frame(to_export_temp)
    to_export=rbind(to_export,to_export_temp)
  }  
  
  
  (P_value = pnorm(Expected_MA, CI_all_chromosomes[,1], Bootstrap_sd))
  print (P_value)
  to_export=cbind(to_export,P_value)
  output<-paste('C:/Users/Waters/Dropbox (MERLAB)/CESRF_2014_adult_sequencing/final_analysis/sliding_window_analysis/',name_output,'_kernel_smoothing_',Nb_bootstrap,'_bootstraps_sigma_',Sigma_sliding_window,'_div',division,'.txt', sep='')
  colnames(to_export)=c("chromosome","position","Fst/Fct","Mean_boostrap","lower_95","upper_95","pvalue")
  print(head(to_export))
  write.table(to_export,file=output,quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
}



#windows()
Dat = read.table("sliding_window_input_2010seg.txt", header=TRUE, sep="\t")
colnames(Dat) = c("Locus","fst", "chromosome", "position")
### Re-order the dataset so that for each chromosome, I have data ordered by the closest to furthest cm from the telomere
Dat = Dat[order(Dat$chromosome, Dat$position),]


par(mfrow=c(6,6), mar=c(2,2,2,2))
plotting_reg_interval(Sigma_sliding_window=3, Nb_bootstrap=1000000, which.chromosome.analysis="all",division=100,name_output="2010SEG")


###############################################################################################################################################################################################

# 2014INT sigma 3 (window=18cM); 1,000,000 reps

par(mfrow=c(6,6), mar=c(2,2,2,2)) ### change to correspond to your data
plotting_reg_interval = function(Sigma_sliding_window=3, Nb_bootstrap=1000000, which.chromosome.analysis="all", which.chromosome.plot="All",division=100,name_output="2014INT") ### these parameters should be specified. If they are not, the default will be used
{
  
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
    plot(Positions_LG, MAs_exp, type="l", ylim=c(min_lim,max_lim), lwd=3, col="red", main = paste(ijk))
    lines(Positions_LG, CI_all[,1], col=1, lwd=2)
    lines(Positions_LG, CI_all[,2], col=rgb(0.3,0.1,0.4, 0.5))
    lines(Positions_LG, CI_all[,3], col=rgb(0.3,0.1,0.4, 0.5))
    if(length(where.NA)!=0) polygon(c(Positions_LG[-where.NA], rev(Positions_LG[-where.NA])), c(CI_all[,2][-where.NA], rev(CI_all[,3][-where.NA])), col=rgb(0.3,0.1,0.4, 0.5))
    if(length(where.NA)==0) polygon(c(Positions_LG, rev(Positions_LG)), c(CI_all[,2], rev(CI_all[,3])), col=rgb(0.3,0.1,0.4, 0.5))
    positions_to_add=Positions_LG
    chromosome=matrix(ijk,length(positions_to_add),1)
    to_export_temp=cbind(chromosome,positions_to_add,MAs_exp,CI_all)
    to_export_temp=as.data.frame(to_export_temp)
    to_export=rbind(to_export,to_export_temp)
  }  
  
  
  (P_value = pnorm(Expected_MA, CI_all_chromosomes[,1], Bootstrap_sd))
  print (P_value)
  to_export=cbind(to_export,P_value)
  output<-paste('C:/Users/Waters/Dropbox (MERLAB)/CESRF_2014_adult_sequencing/final_analysis/sliding_window_analysis/',name_output,'_kernel_smoothing_',Nb_bootstrap,'_bootstraps_sigma_',Sigma_sliding_window,'_div',division,'.txt', sep='')
  colnames(to_export)=c("chromosome","position","Fst/Fct","Mean_boostrap","lower_95","upper_95","pvalue")
  print(head(to_export))
  write.table(to_export,file=output,quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
}



#windows()
Dat = read.table("sliding_window_input_2014int.txt", header=TRUE, sep="\t")
colnames(Dat) = c("Locus","fst", "chromosome", "position")
### Re-order the dataset so that for each chromosome, I have data ordered by the closest to furthest cm from the telomere
Dat = Dat[order(Dat$chromosome, Dat$position),]


par(mfrow=c(6,6), mar=c(2,2,2,2))
plotting_reg_interval(Sigma_sliding_window=3, Nb_bootstrap=1000000, which.chromosome.analysis="all",division=100,name_output="2014INT")


###############################################################################################################################################################################################

# 2014SEG sigma 3 (window=18cM); 1,000,000 reps

par(mfrow=c(6,6), mar=c(2,2,2,2)) ### change to correspond to your data
plotting_reg_interval = function(Sigma_sliding_window=3, Nb_bootstrap=1000000, which.chromosome.analysis="all", which.chromosome.plot="All",division=100,name_output="2014SEG") ### these parameters should be specified. If they are not, the default will be used
{
  
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
    plot(Positions_LG, MAs_exp, type="l", ylim=c(min_lim,max_lim), lwd=3, col="red", main = paste(ijk))
    lines(Positions_LG, CI_all[,1], col=1, lwd=2)
    lines(Positions_LG, CI_all[,2], col=rgb(0.3,0.1,0.4, 0.5))
    lines(Positions_LG, CI_all[,3], col=rgb(0.3,0.1,0.4, 0.5))
    if(length(where.NA)!=0) polygon(c(Positions_LG[-where.NA], rev(Positions_LG[-where.NA])), c(CI_all[,2][-where.NA], rev(CI_all[,3][-where.NA])), col=rgb(0.3,0.1,0.4, 0.5))
    if(length(where.NA)==0) polygon(c(Positions_LG, rev(Positions_LG)), c(CI_all[,2], rev(CI_all[,3])), col=rgb(0.3,0.1,0.4, 0.5))
    positions_to_add=Positions_LG
    chromosome=matrix(ijk,length(positions_to_add),1)
    to_export_temp=cbind(chromosome,positions_to_add,MAs_exp,CI_all)
    to_export_temp=as.data.frame(to_export_temp)
    to_export=rbind(to_export,to_export_temp)
  }  
  
  
  (P_value = pnorm(Expected_MA, CI_all_chromosomes[,1], Bootstrap_sd))
  print (P_value)
  to_export=cbind(to_export,P_value)
  output<-paste('C:/Users/Waters/Dropbox (MERLAB)/CESRF_2014_adult_sequencing/final_analysis/sliding_window_analysis/',name_output,'_kernel_smoothing_',Nb_bootstrap,'_bootstraps_sigma_',Sigma_sliding_window,'_div',division,'.txt', sep='')
  colnames(to_export)=c("chromosome","position","Fst/Fct","Mean_boostrap","lower_95","upper_95","pvalue")
  print(head(to_export))
  write.table(to_export,file=output,quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
}



#windows()
Dat = read.table("sliding_window_input_2014seg.txt", header=TRUE, sep="\t")
colnames(Dat) = c("Locus","fst", "chromosome", "position")
### Re-order the dataset so that for each chromosome, I have data ordered by the closest to furthest cm from the telomere
Dat = Dat[order(Dat$chromosome, Dat$position),]


par(mfrow=c(6,6), mar=c(2,2,2,2))
plotting_reg_interval(Sigma_sliding_window=3, Nb_bootstrap=1000000, which.chromosome.analysis="all",division=100,name_output="2014SEG")
