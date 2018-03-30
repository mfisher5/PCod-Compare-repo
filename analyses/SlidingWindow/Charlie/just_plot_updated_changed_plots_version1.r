###Calculate the number of data points (loci) within each window for each chromosome


marker_data = read.table("C:/Users/Waters/Dropbox (MERLAB)/CESRF_2014_adult_sequencing/final_analysis/sliding_window_analysis/test_density_average.txt", header=TRUE, sep="\t")

head(marker_data)
colnames(marker_data) = c("Locus","fst", "chromosome", "position")
### Re-order the dataset so that for each chromosome, I have data ordered by the closest to furthest cm from the telomere
marker_data = marker_data[order(marker_data$chromosome, marker_data$position),]

window_size=3
divisions=100
Num_marker_data = nrow(marker_data)
Num_chromosome = unique(marker_data$chromosome)
chromosomes_all= rep(Num_chromosome,each=100)
Positions_all=c()
markers_per_window=c()
par(mfrow=c(1,1), mar=c(2,2,2,2))

for (ijx in Num_chromosome){
  
  Data.chrom = subset(marker_data, marker_data$chromosome==ijx)
  Size_windows = 3*window_size  # sliding window is 3 times the value of sigma in each direction
  Beg_position_cM=min(Data.chrom$position)+Size_windows  ## calculations cannot start at 0, because it would be biased no values between -6 and 0cM for example
  End_position_cM=max(Data.chrom$position)-Size_windows
  Beg_position = which(Data.chrom$position >=Beg_position_cM)[1] ### start at the first position as calculated before or at the first position where there is a data point (as sometimes, a large portion of the chromosome is made of duplicated loci --> not kept for population studies)
  End_position = which(Data.chrom$position >= (max(End_position_cM)))[1]
  Positions_LG=seq(from=Beg_position_cM,to=End_position_cM,length.out=divisions)  ## these position will be equidistant
  Positions_all=c(Positions_all,Positions_LG)
  markers_each_chrom=c()
  for (x in 1:length(Positions_LG)){
    MA_window = which(Data.chrom$position >= (Positions_LG[x]-Size_windows) & Data.chrom$position <= (Positions_LG[x]+Size_windows))
    num_markers = length(MA_window)
    #markers_each_chrom=c(markers_each_chrom,num_markers)
    markers_per_window = c(markers_per_window,num_markers)
    
  }
  #plot(Positions_LG,markers_each_chrom,main=ijx)
  #lines(Positions_LG,markers_each_chrom,type='l')
}

output=c()   ##be sure to clear output from previous runs
output=cbind(chromosomes_all,Positions_all,markers_per_window)
colnames(output)=c("chrom","position","num_markers")
markers_per_window




###Now these are Marine's (and Kot's) functions to plot the moving average of Fst
par(mar=c(6,6,3,3))

#### 2002INT
Data_to_plot_2002INT<-read.table('C:/Users/Waters/Dropbox (MERLAB)/Adult_traits_paper/Figures/Sliding_window/2002INT_kernel_smoothing_1e+06_bootstraps_sigma_3_div100.txt', header=TRUE, sep="\t")

Data_to_plot_2002INT <- cbind(Data_to_plot_2002INT, markers_per_window)

just_plot_2002INT = function(Nb_bootstrap=1000000, which.chromosome.analysis="all", which.chromosome.plot="all",name="1998_2002INT")
{
  
  if(which.chromosome.analysis!="all") Data.analysis_2002INT = subset(Data_to_plot_2002INT, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis=="all") Data.analysis_2002INT = Data_to_plot_2002INT
  
  Nb_chromosome = unique(Data.analysis_2002INT$chromosome)
  for (ijk in Nb_chromosome)
  {
    
    Data.part = subset(Data_to_plot_2002INT, subset=chromosome %in% ijk)
    print(head(Data.part))
    min_lim=min(Data.part[,3:6],na.rm=TRUE)-0.02
    max_lim=max(Data.part[,3:6],na.rm=TRUE)+0.1
    where.NA = which(is.na(Data.part[,3]))
    plot(Data.part[,2], Data.part[,3],type="l",ylim=c(-0.02,0.2), xlim=c(0,150), lwd=3, col="red", ylab='',xlab='',cex.lab=2, cex.axis=2)
    title(expression(paste("F"[1], "Wild")),cex.main=2.5,line=-1)
    lines(Data.part[,2], Data.part[,5], col=rgb(0.193,0.205,0.205,0.25))
    lines(Data.part[,2], Data.part[,6], col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)!=0) polygon(c(Data.part[,2][-where.NA], rev(Data.part[,2][-where.NA])), c(Data.part[,5][-where.NA], rev(Data.part[,6][-where.NA])), col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)==0) polygon(c(Data.part[,2], rev(Data.part[,2])), c(Data.part[,5], rev(Data.part[,6])), col=rgb(0.193,0.205,0.205,0.25))    
    abline(h=0)
    #par(new=T)
    #plot(Data.part[,2],Data.part[,8],ylim=c(-10,1.1*max(Data.part[,8])),xlim=c(0,155),col="purple",type='l',lwd=2,axes=F,ylab='',xlab='')
    #axis(4,at=c(0,10,20,30,40,50,60,70),labels=TRUE,cex.axis=1.5,cex=1.5)
    #mtext("Loci per Window",side=4,line=2.5,cex=1.5)
  
  }  
  
}

#head(Data_to_plot_2002INT)
#par(mfrow=c(6,6), mar=c(2,2,2,2))
#just_plot_2002INT(Nb_bootstrap=1000000, which.chromosome.analysis="all",which.chromosome.plot="all" ,name="1998_2002INT")


#### 2002SEG
Data_to_plot_2002SEG<-read.table('C:/Users/Waters/Dropbox (MERLAB)/Adult_traits_paper/Figures/Sliding_window/2002SEG_kernel_smoothing_1e+06_bootstraps_sigma_3_div100.txt', header=TRUE, sep="\t")
Data_to_plot_2002SEG <- cbind(Data_to_plot_2002SEG, markers_per_window)

just_plot_2002SEG = function(Nb_bootstrap=1000000, which.chromosome.analysis="all", which.chromosome.plot="all",name="1998_2002SEG")
{
  
  if(which.chromosome.analysis!="all") Data.analysis_2002SEG = subset(Data_to_plot_2002SEG, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis=="all") Data.analysis_2002SEG = Data_to_plot_2002SEG
  
  Nb_chromosome = unique(Data.analysis_2002SEG$chromosome)
  for (ijk in Nb_chromosome)
  {
    
    Data.part = subset(Data_to_plot_2002SEG, subset=chromosome %in% ijk)
    print(head(Data.part))
    min_lim=min(Data.part[,3:6],na.rm=TRUE)-0.02
    max_lim=max(Data.part[,3:6],na.rm=TRUE)+0.1
    where.NA = which(is.na(Data.part[,3]))
    plot(Data.part[,2], Data.part[,3],type="l",ylim=c(-0.02,0.2), xlim=c(0,150), lwd=3, col="red",ylab='', xlab='',cex.lab=2, cex.axis=2)
    title(expression(paste("F"[1], "Hatchery")),cex.main=2.5,line=-1)
    lines(Data.part[,2], Data.part[,5], col=rgb(0.193,0.205,0.205,0.25))
    lines(Data.part[,2], Data.part[,6], col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)!=0) polygon(c(Data.part[,2][-where.NA], rev(Data.part[,2][-where.NA])), c(Data.part[,5][-where.NA], rev(Data.part[,6][-where.NA])), col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)==0) polygon(c(Data.part[,2], rev(Data.part[,2])), c(Data.part[,5], rev(Data.part[,6])), col=rgb(0.193,0.205,0.205,0.25))    
    abline(h=0)
    #par(new=T)
    #plot(Data.part[,2],Data.part[,8],ylim=c(-10,1.1*max(Data.part[,8])),xlim=c(0,155),col="purple",type='l',lwd=2,axes=F,ylab='',xlab='')
    #axis(4,at=c(0,10,20,30,40,50,60,70),labels=TRUE,cex.axis=1.5,cex=1.5)
    #mtext("Loci per Window",side=4,line=2.5,cex=1.5)
    
  }  
  
}

#head(Data_to_plot_2002SEG)
#par(mfrow=c(6,6), mar=c(2,2,2,2))
#just_plot_2002SEG(Nb_bootstrap=1000000, which.chromosome.analysis="all",which.chromosome.plot="all" ,name="1998_2002SEG")


#### 2006INT
Data_to_plot_2006INT<-read.table('C:/Users/Waters/Dropbox (MERLAB)/Adult_traits_paper/Figures/Sliding_window/2006INT_kernel_smoothing_1e+06_bootstraps_sigma_3_div100.txt', header=TRUE, sep="\t")
Data_to_plot_2006INT <- cbind(Data_to_plot_2006INT, markers_per_window)

just_plot_2006INT = function(Nb_bootstrap=1000000, which.chromosome.analysis="all", which.chromosome.plot="all",name="1998_2006INT")
{
  
  if(which.chromosome.analysis!="all") Data.analysis_2006INT = subset(Data_to_plot_2006INT, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis=="all") Data.analysis_2006INT = Data_to_plot_2006INT
  
  Nb_chromosome = unique(Data.analysis_2006INT$chromosome)
  for (ijk in Nb_chromosome)
  {
    
    Data.part = subset(Data_to_plot_2006INT, subset=chromosome %in% ijk)
    print(head(Data.part))
    min_lim=min(Data.part[,3:6],na.rm=TRUE)-0.02
    max_lim=max(Data.part[,3:6],na.rm=TRUE)+0.1
    where.NA = which(is.na(Data.part[,3]))
    plot(Data.part[,2], Data.part[,3],type="l",ylim=c(-0.02,0.2), xlim=c(0,150), lwd=3, col="red", ylab='',xlab='',cex.lab=2, cex.axis=2)
    title(expression(paste("F"[2], "INT")),cex.main=2.5,line=-1)
    lines(Data.part[,2], Data.part[,5], col=rgb(0.193,0.205,0.205,0.25))
    lines(Data.part[,2], Data.part[,6], col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)!=0) polygon(c(Data.part[,2][-where.NA], rev(Data.part[,2][-where.NA])), c(Data.part[,5][-where.NA], rev(Data.part[,6][-where.NA])), col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)==0) polygon(c(Data.part[,2], rev(Data.part[,2])), c(Data.part[,5], rev(Data.part[,6])), col=rgb(0.193,0.205,0.205,0.25))    
    abline(h=0)
    #par(new=T)
    #plot(Data.part[,2],Data.part[,8],ylim=c(-10,1.1*max(Data.part[,8])),xlim=c(0,155),col="purple",type='l',lwd=2,axes=F,ylab='',xlab='')
    #axis(4,at=c(0,10,20,30,40,50,60,70),labels=TRUE,cex.axis=1.5,cex=1.5)
    #mtext("Loci per Window",side=4,line=2.5,cex=1.5)
    
  }  
  
}

#head(Data_to_plot_2006INT)
#par(mfrow=c(6,6), mar=c(2,2,2,2))
#just_plot_2006INT(Nb_bootstrap=1000000, which.chromosome.analysis="all",which.chromosome.plot="all" ,name="1998_2006INT")


#### 2006SEG
Data_to_plot_2006SEG<-read.table('C:/Users/Waters/Dropbox (MERLAB)/Adult_traits_paper/Figures/Sliding_window/2006SEG_kernel_smoothing_1e+06_bootstraps_sigma_3_div100.txt', header=TRUE, sep="\t")
Data_to_plot_2006SEG <- cbind(Data_to_plot_2006SEG, markers_per_window)

just_plot_2006SEG = function(Nb_bootstrap=1000000, which.chromosome.analysis="all", which.chromosome.plot="all",name="1998_2006SEG")
{
  
  if(which.chromosome.analysis!="all") Data.analysis_2006SEG = subset(Data_to_plot_2006SEG, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis=="all") Data.analysis_2006SEG = Data_to_plot_2006SEG
  
  Nb_chromosome = unique(Data.analysis_2006SEG$chromosome)
  for (ijk in Nb_chromosome)
  {
    
    Data.part = subset(Data_to_plot_2006SEG, subset=chromosome %in% ijk)
    print(head(Data.part))
    min_lim=min(Data.part[,3:6],na.rm=TRUE)-0.02
    max_lim=max(Data.part[,3:6],na.rm=TRUE)+0.1
    where.NA = which(is.na(Data.part[,3]))
    plot(Data.part[,2], Data.part[,3],type="l",ylim=c(-0.02,0.2), xlim=c(0,150), lwd=3, col="red", ylab='', xlab='',cex.lab=2, cex.axis=2)
    title(expression(paste("F"[2], "SEG")),cex.main=2.5,line=-1)
    lines(Data.part[,2], Data.part[,5], col=rgb(0.193,0.205,0.205,0.25))
    lines(Data.part[,2], Data.part[,6], col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)!=0) polygon(c(Data.part[,2][-where.NA], rev(Data.part[,2][-where.NA])), c(Data.part[,5][-where.NA], rev(Data.part[,6][-where.NA])), col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)==0) polygon(c(Data.part[,2], rev(Data.part[,2])), c(Data.part[,5], rev(Data.part[,6])), col=rgb(0.193,0.205,0.205,0.25))    
    abline(h=0)
    #par(new=T)
    #plot(Data.part[,2],Data.part[,8],ylim=c(-10,1.1*max(Data.part[,8])),xlim=c(0,155),col="purple",type='l',lwd=2,axes=F,ylab='',xlab='')
    #axis(4,at=c(0,10,20,30,40,50,60,70),labels=TRUE,cex.axis=1.5,cex=1.5)
    #mtext("Loci per Window",side=4,line=2.5,cex=1.5)
    
  }  
  
}

#head(Data_to_plot_2006SEG)
#par(mfrow=c(6,6), mar=c(2,2,2,2))
#just_plot_2006SEG(Nb_bootstrap=1000000, which.chromosome.analysis="all",which.chromosome.plot="all" ,name="1998_2006SEG")


#### 2010INT
Data_to_plot_2010INT<-read.table('C:/Users/Waters/Dropbox (MERLAB)/Adult_traits_paper/Figures/Sliding_window/2010INT_kernel_smoothing_1e+06_bootstraps_sigma_3_div100.txt', header=TRUE, sep="\t")
Data_to_plot_2010INT <- cbind(Data_to_plot_2010INT, markers_per_window)

just_plot_2010INT = function(Nb_bootstrap=1000000, which.chromosome.analysis="all", which.chromosome.plot="all",name="1998_2010INT")
{
  
  if(which.chromosome.analysis!="all") Data.analysis_2010INT = subset(Data_to_plot_2010INT, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis=="all") Data.analysis_2010INT = Data_to_plot_2010INT
  
  Nb_chromosome = unique(Data.analysis_2010INT$chromosome)
  for (ijk in Nb_chromosome)
  {
    
    Data.part = subset(Data_to_plot_2010INT, subset=chromosome %in% ijk)
    print(head(Data.part))
    min_lim=min(Data.part[,3:6],na.rm=TRUE)-0.02
    max_lim=max(Data.part[,3:6],na.rm=TRUE)+0.1
    where.NA = which(is.na(Data.part[,3]))
    plot(Data.part[,2], Data.part[,3],type="l",ylim=c(-0.02,0.2), xlim=c(0,150), lwd=3, col="red", ylab='', xlab='',cex.lab=2, cex.axis=2)
    title(expression(paste("F"[3], "INT")),cex.main=2.5,line=-1)
    lines(Data.part[,2], Data.part[,5], col=rgb(0.193,0.205,0.205,0.25))
    lines(Data.part[,2], Data.part[,6], col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)!=0) polygon(c(Data.part[,2][-where.NA], rev(Data.part[,2][-where.NA])), c(Data.part[,5][-where.NA], rev(Data.part[,6][-where.NA])), col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)==0) polygon(c(Data.part[,2], rev(Data.part[,2])), c(Data.part[,5], rev(Data.part[,6])), col=rgb(0.193,0.205,0.205,0.25))    
    abline(h=0)
    #par(new=T)
    #plot(Data.part[,2],Data.part[,8],ylim=c(-10,1.1*max(Data.part[,8])),xlim=c(0,155),col="purple",type='l',lwd=2,axes=F,ylab='',xlab='')
    #axis(4,at=c(0,10,20,30,40,50,60,70),labels=TRUE,cex.axis=1.5,cex=1.5)
    #mtext("Loci per Window",side=4,line=2.5,cex=1.5)
    
  }  
  
}

#head(Data_to_plot_2010INT)
#par(mfrow=c(6,6), mar=c(2,2,2,2))
#just_plot_2010INT(Nb_bootstrap=1000000, which.chromosome.analysis="all",which.chromosome.plot="all" ,name="1998_2010INT")


#### 2010SEG
Data_to_plot_2010SEG<-read.table('C:/Users/Waters/Dropbox (MERLAB)/Adult_traits_paper/Figures/Sliding_window/2010SEG_kernel_smoothing_1e+06_bootstraps_sigma_3_div100.txt', header=TRUE, sep="\t")
Data_to_plot_2010SEG <- cbind(Data_to_plot_2010SEG, markers_per_window)

just_plot_2010SEG = function(Nb_bootstrap=1000000, which.chromosome.analysis="all", which.chromosome.plot="all",name="1998_2010SEG")
{
  
  if(which.chromosome.analysis!="all") Data.analysis_2010SEG = subset(Data_to_plot_2010SEG, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis=="all") Data.analysis_2010SEG = Data_to_plot_2010SEG
  
  Nb_chromosome = unique(Data.analysis_2010SEG$chromosome)
  for (ijk in Nb_chromosome)
  {
    
    Data.part = subset(Data_to_plot_2010SEG, subset=chromosome %in% ijk)
    print(head(Data.part))
    min_lim=min(Data.part[,3:6],na.rm=TRUE)-0.02
    max_lim=max(Data.part[,3:6],na.rm=TRUE)+0.1
    where.NA = which(is.na(Data.part[,3]))
    plot(Data.part[,2], Data.part[,3],type="l",ylim=c(-0.02,0.2), xlim=c(0,150), lwd=3, col="red", ylab='', xlab='',cex.lab=2, cex.axis=2)
    title(expression(paste("F"[3], "SEG")),cex.main=2.5,line=-1)
    lines(Data.part[,2], Data.part[,5], col=rgb(0.193,0.205,0.205,0.25))
    lines(Data.part[,2], Data.part[,6], col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)!=0) polygon(c(Data.part[,2][-where.NA], rev(Data.part[,2][-where.NA])), c(Data.part[,5][-where.NA], rev(Data.part[,6][-where.NA])), col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)==0) polygon(c(Data.part[,2], rev(Data.part[,2])), c(Data.part[,5], rev(Data.part[,6])), col=rgb(0.193,0.205,0.205,0.25))    
    abline(h=0)
    #par(new=T)
    #plot(Data.part[,2],Data.part[,8],ylim=c(-10,1.1*max(Data.part[,8])),xlim=c(0,155),col="purple",type='l',lwd=2,axes=F,ylab='',xlab='')
    #axis(4,at=c(0,10,20,30,40,50,60,70),labels=TRUE,cex.axis=1.5,cex=1.5)
    #mtext("Loci per Window",side=4,line=2.5,cex=1.5)
    
  }  
  
}


#head(Data_to_plot_2010SEG)
#par(mfrow=c(6,6), mar=c(2,2,2,2))
#just_plot_2010SEG(Nb_bootstrap=1000000, which.chromosome.analysis="all",which.chromosome.plot="all" ,name="1998_2010SEG")


#### 2014INT
Data_to_plot_2014INT<-read.table('C:/Users/Waters/Dropbox (MERLAB)/Adult_traits_paper/Figures/Sliding_window/2014INT_kernel_smoothing_1e+06_bootstraps_sigma_3_div100.txt', header=TRUE, sep="\t")
Data_to_plot_2014INT <- cbind(Data_to_plot_2014INT, markers_per_window)

just_plot_2014INT = function(Nb_bootstrap=1000000, which.chromosome.analysis="all", which.chromosome.plot="all",name="1998_2014INT")
{
  
  if(which.chromosome.analysis!="all") Data.analysis_2014INT = subset(Data_to_plot_2014INT, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis=="all") Data.analysis_2014INT = Data_to_plot_2014INT
  
  Nb_chromosome = unique(Data.analysis_2014INT$chromosome)
  for (ijk in Nb_chromosome)
  {
    
    Data.part = subset(Data_to_plot_2014INT, subset=chromosome %in% ijk)
    print(head(Data.part))
    min_lim=min(Data.part[,3:6],na.rm=TRUE)-0.02
    max_lim=max(Data.part[,3:6],na.rm=TRUE)+0.1
    where.NA = which(is.na(Data.part[,3]))
    plot(Data.part[,2], Data.part[,3],type="l",ylim=c(-0.02,0.2), xlim=c(0,150), lwd=3, col="red", ylab='', xlab='',cex.lab=2, cex.axis=2)
    title(expression(paste("F"[4], "INT")),cex.main=2.5,line=-1)
    lines(Data.part[,2], Data.part[,5], col=rgb(0.193,0.205,0.205,0.25))
    lines(Data.part[,2], Data.part[,6], col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)!=0) polygon(c(Data.part[,2][-where.NA], rev(Data.part[,2][-where.NA])), c(Data.part[,5][-where.NA], rev(Data.part[,6][-where.NA])), col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)==0) polygon(c(Data.part[,2], rev(Data.part[,2])), c(Data.part[,5], rev(Data.part[,6])), col=rgb(0.193,0.205,0.205,0.25))    
    abline(h=0)
    #par(new=T)
    #plot(Data.part[,2],Data.part[,8],ylim=c(-10,1.1*max(Data.part[,8])),xlim=c(0,155),col="purple",type='l',lwd=2,axes=F,ylab='',xlab='')
    #axis(4,at=c(0,10,20,30,40,50,60,70),labels=TRUE,cex.axis=1.5,cex=1.5)
    #mtext("Loci per Window",side=4,line=2.5,cex=1.5)
    
  }  
  
}

#head(Data_to_plot_2014INT)
#par(mfrow=c(6,6), mar=c(2,2,2,2))
#just_plot_2014INT(Nb_bootstrap=1000000, which.chromosome.analysis="all",which.chromosome.plot="all" ,name="1998_2014INT")


#### 2014SEG
Data_to_plot_2014SEG<-read.table('C:/Users/Waters/Dropbox (MERLAB)/Adult_traits_paper/Figures/Sliding_window/2014SEG_kernel_smoothing_1e+06_bootstraps_sigma_3_div100.txt', header=TRUE, sep="\t")
Data_to_plot_2014SEG <- cbind(Data_to_plot_2014SEG, markers_per_window)

just_plot_2014SEG = function(Nb_bootstrap=1000000, which.chromosome.analysis="all", which.chromosome.plot="all",name="1998_2014SEG")
{
  
  if(which.chromosome.analysis!="all") Data.analysis_2014SEG = subset(Data_to_plot_2014SEG, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis=="all") Data.analysis_2014SEG = Data_to_plot_2014SEG
  
  Nb_chromosome = unique(Data.analysis_2014SEG$chromosome)
  for (ijk in Nb_chromosome)
  {
    
    Data.part = subset(Data_to_plot_2014SEG, subset=chromosome %in% ijk)
    print(head(Data.part))
    min_lim=min(Data.part[,3:6],na.rm=TRUE)-0.02
    max_lim=max(Data.part[,3:6],na.rm=TRUE)+0.1
    where.NA = which(is.na(Data.part[,3]))
    plot(Data.part[,2], Data.part[,3],type="l",ylim=c(-0.02,0.2), xlim=c(0,150), lwd=3, col="red", ylab='', xlab='',cex.lab=2, cex.axis=2)
    title(expression(paste("F"[4], "SEG")),cex.main=2.5,line=-1)
    lines(Data.part[,2], Data.part[,5], col=rgb(0.193,0.205,0.205,0.25))
    lines(Data.part[,2], Data.part[,6], col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)!=0) polygon(c(Data.part[,2][-where.NA], rev(Data.part[,2][-where.NA])), c(Data.part[,5][-where.NA], rev(Data.part[,6][-where.NA])), col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)==0) polygon(c(Data.part[,2], rev(Data.part[,2])), c(Data.part[,5], rev(Data.part[,6])), col=rgb(0.193,0.205,0.205,0.25))    
    abline(h=0)
    #par(new=T)
    #plot(Data.part[,2],Data.part[,8],ylim=c(-10,1.1*max(Data.part[,8])),xlim=c(0,155),col="purple",type='l',lwd=2,axes=F,ylab='',xlab='')
    #axis(4,at=c(0,10,20,30,40,50,60,70),labels=TRUE,cex.axis=1.5,cex=1.5)
    #mtext("Loci per Window",side=4,line=2.5,cex=1.5)
    
  }  
  
}


#head(Data_to_plot_2014SEG)
#par(mfrow=c(6,6), mar=c(2,2,2,2))
#just_plot_2014SEG(Nb_bootstrap=1000000, which.chromosome.analysis="all",which.chromosome.plot="all" ,name="1998_2014SEG")





par()  ##view default par settings
opar<-par()   ##save a copy of default par settings so that I can revert back to them while plotting
par<-opar

## Set plot window specifications
par(mfrow=c(2,4), mar=c(6,5,5,4),oma=c(5,1,3,1),mai=c(0.5,1,0.5,0.15))

## Read in data
data_1998_2002INT_Fst <- read.table('C:/Users/Waters/Dropbox (MERLAB)/Adult_traits_paper/Figures/Sliding_window/sliding_window_input_2002int.txt', header=TRUE, sep="\t")
data_1998_2006INT_Fst <- read.table('C:/Users/Waters/Dropbox (MERLAB)/Adult_traits_paper/Figures/Sliding_window/sliding_window_input_2006int.txt', header=TRUE, sep="\t")
data_1998_2010INT_Fst <- read.table('C:/Users/Waters/Dropbox (MERLAB)/Adult_traits_paper/Figures/Sliding_window/sliding_window_input_2010int.txt', header=TRUE, sep="\t")
data_1998_2014INT_Fst <- read.table('C:/Users/Waters/Dropbox (MERLAB)/Adult_traits_paper/Figures/Sliding_window/sliding_window_input_2014int.txt', header=TRUE, sep="\t")
data_1998_2002SEG_Fst <- read.table('C:/Users/Waters/Dropbox (MERLAB)/Adult_traits_paper/Figures/Sliding_window/sliding_window_input_2002seg.txt', header=TRUE, sep="\t")
data_1998_2006SEG_Fst <- read.table('C:/Users/Waters/Dropbox (MERLAB)/Adult_traits_paper/Figures/Sliding_window/sliding_window_input_2006seg.txt', header=TRUE, sep="\t")
data_1998_2010SEG_Fst <- read.table('C:/Users/Waters/Dropbox (MERLAB)/Adult_traits_paper/Figures/Sliding_window/sliding_window_input_2010seg.txt', header=TRUE, sep="\t")
data_1998_2014SEG_Fst <- read.table('C:/Users/Waters/Dropbox (MERLAB)/Adult_traits_paper/Figures/Sliding_window/sliding_window_input_2014seg.txt', header=TRUE, sep="\t")


#Fst values for all 4214 mapped loci
#Fst_all_loci <- read.csv('C:/Users/Waters/Dropbox (MERLAB)/Adult_traits_paper/Figures/Sliding_window/4gens_Fst_4214mapped.csv', header=TRUE)


###############
all_outliers<-read.csv('C:/Users/Waters/Dropbox (MERLAB)/Adult_traits_paper/Figures/Sliding_window/All_final_outliers_mapped.csv', header=TRUE)

bayescan_outliers <- all_outliers[all_outliers[,"Method"]=="Bayescan",]

ftemp_outliers_SEG <- all_outliers[all_outliers[,"Method"]=="FtempSEG",]

ftemp_outliers_INT <- all_outliers[all_outliers[,"Method"]=="FtempINT",]

#dapc1_structure_loci <- all_outliers[all_outliers[,"Method"]=="DAPC",]

#bayescan_FtempSEG_overlap<-read.csv('C:/Users/Waters/Dropbox (MERLAB)/CESRF_2014_adult_sequencing/final_analysis/sliding_window_analysis/Bayescan_FtempSEG_overlap.csv', header=TRUE)

#bayescan_DAPC_ftempSEG<-read.csv('C:/Users/Waters/Dropbox (MERLAB)/Analysis_locus_filter/sliding_window_analysis/Bayescan_DAPC_FtempSEG_overlap.csv', header=TRUE)
####################

for (ijx in Num_chromosome[1:34]) {
  ijx<-"Ots12"
  par(mfrow=c(2,4),mar=c(6,6,5,6),oma=c(8,1,3,4),mai=c(0.25,0.35,0.15,0.35))
  
  #par(mfrow=c(2,1),mar=c(6,6,5,4),oma=c(4,1,1,1),mai=c(0.25,0.5,0.25,0.25))
  
  data_2002INT <- data_1998_2002INT_Fst[data_1998_2002INT_Fst[,"chromosome"]==ijx,]
  chrom_2002INT <- just_plot_2002INT(Nb_bootstrap=1000000, which.chromosome.analysis=ijx,which.chromosome.plot=ijx ,name="1998_2002INT")
  chrom_loci <- Fst_all_loci[Fst_all_loci[,"Chrom"]==ijx,]
  #points(chrom_loci$Position,chrom_loci$fst2002int,type="p",col="black",pch=16,cex=1.2)
  chrom_bayes <-bayescan_outliers[bayescan_outliers[,"Chrom"]==ijx,]
  points(chrom_bayes$position,chrom_bayes$fst2002int,type='p',col="blue",pch=15,cex=3)
  chrom_dapc <- dapc1_structure_loci[dapc1_structure_loci[,"Chrom"]==ijx,]
  points(chrom_dapc$position,chrom_dapc$fst2002int,type='p',col="darkgreen",pch=16,cex=2.5)
  chrom_ftempINT <- ftemp_outliers_INT[ftemp_outliers_INT[,"Chrom"]==ijx,]
  points(chrom_ftempINT$position,chrom_ftempINT$fst2002int,type='p',col="orange",pch=17,cex=2.5)
  #chrom_bayes_dapc <- bayescan_DAPC_ftempSEG[bayescan_DAPC_ftempSEG[,"Chrom"]==ijx,]
  #points(chrom_bayes_dapc$position,chrom_bayes_dapc$fst2002int,type='p',col="orange",pch=18,cex=3)
  
  #add predictor locus for return time that Marine found on Ots12
  #points(20.7,-0.0138, type='p',col="black",pch=18,cex=4)
  
  
  #par(new=T)
  #Data.part_2002INT = subset(Data_to_plot_2002INT, subset=chromosome==ijx)
  #plot(Data.part_2002INT[,2],Data.part_2002INT[,8],ylim=c(-10,100),xlim=c(0,150),xaxs='i',col="purple",type='l',lwd=2,axes=F,ylab='',xlab='')
  #axis(4,at=c(0,10,20,30,40,50,60,70,80,90,100),labels=FALSE,cex.axis=1.5,cex=1.5)
  #axis(4,at=c(0,20,40,60,80,100),labels=TRUE,cex.axis=2,cex=1.5)
  
############################################################################################################

  data_2006INT <- data_1998_2006INT_Fst[data_1998_2006INT_Fst[,"chromosome"]==ijx,]
  chrom_2006INT <- just_plot_2006INT(Nb_bootstrap=1000000, which.chromosome.analysis=ijx ,which.chromosome.plot=ijx ,name="1998_2006INT")
  #points(chrom_loci$Position,chrom_loci$fst2006int,type="p",col="black",pch=16,cex=1.2)
  points(chrom_bayes$position,chrom_bayes$fst2006int,type='p',col="blue",pch=15,cex=3)
  points(chrom_dapc$position,chrom_dapc$fst2006int,type='p',col="darkgreen",pch=16,cex=2.5)
  points(chrom_ftempINT$position,chrom_ftempINT$fst2006int,type='p',col="orange",pch=17,cex=2.5)
  #points(chrom_bayes_dapc$position,chrom_bayes_dapc$fst2006int,type='p',col="orange",pch=18,cex=3)

  #add predictor locus for return time that Marine found on Ots12
  #points(20.7,0.012, type='p',col="black",pch=18,cex=4)
  
  
  #par(new=T)
  #Data.part_2006INT = subset(Data_to_plot_2006INT, subset=chromosome==ijx)
  #plot(Data.part_2006INT[,2],Data.part_2006INT[,8],ylim=c(-10,100),xlim=c(0,150),xaxs='i',col="purple",type='l',lwd=2,axes=F,ylab='',xlab='')
  #axis(4,at=c(0,10,20,30,40,50,60,70,80,90,100),labels=FALSE,cex.axis=1.5,cex=1.5)
  #axis(4,at=c(0,20,40,60,80,100),labels=TRUE,cex.axis=2,cex=1.5)
  
############################################################################################################

  data_2010INT <- data_1998_2010INT_Fst[data_1998_2010INT_Fst[,"chromosome"]==ijx,]
  chrom_2010INT <- just_plot_2010INT(Nb_bootstrap=1000000, which.chromosome.analysis=ijx ,which.chromosome.plot=ijx ,name="1998_2010INT")
  #points(chrom_loci$Position,chrom_loci$fst2010int,type="p",col="black",pch=16,cex=1.2)
  points(chrom_bayes$position,chrom_bayes$fst2010int,type='p',col="blue",pch=15,cex=3)
  points(chrom_dapc$position,chrom_dapc$fst2010int,type='p',col="darkgreen",pch=16,cex=2.5)
  points(chrom_ftempINT$position,chrom_ftempINT$fst2010int,type='p',col="orange",pch=17,cex=2.5)
  #points(chrom_bayes_dapc$position,chrom_bayes_dapc$fst2010int,type='p',col="orange",pch=18,cex=3)

  #add predictor locus for return time that Marine found on Ots12
  #points(20.7,-0.0091, type='p',col="black",pch=18,cex=4)
  
  par(new=T)
  Data.part_2010INT = subset(Data_to_plot_2010INT, subset=chromosome==ijx)
  plot(Data.part_2010INT[,2],Data.part_2010INT[,8],ylim=c(-10,100),xlim=c(0,150),xaxs='i',col="purple",type='l',lwd=2,axes=F,ylab='',xlab='')
  #axis(4,at=c(0,10,20,30,40,50,60,70,80,90,100),labels=FALSE,cex.axis=1.5,cex=1.5)
  #axis(4,at=c(0,20,40,60,80,100),labels=TRUE,cex.axis=2,cex=1.5)
  #legend("topright",legend=c("Temporal Outliers","Bayescan Outliers","DAPC Structure Loci","Moving Average", "# Loci per Window"),pch=c(17,15,16,NA,NA),cex=1.5,pt.cex=2,lwd=2,lty=c(NA,NA,NA,1,1),col=c("green","blue","navy","red","purple"))
  
############################################################################################################

  data_2014INT <- data_1998_2014INT_Fst[data_1998_2014INT_Fst[,"chromosome"]==ijx,]
  chrom_2014INT <- just_plot_2014INT(Nb_bootstrap=1000000, which.chromosome.analysis=ijx ,which.chromosome.plot=ijx ,name="1998_2014INT")
  #points(chrom_loci$Position,chrom_loci$fst2014int,type="p",col="black",pch=16,cex=1.2)
  points(chrom_bayes$position,chrom_bayes$fst2014int,type='p',col="blue",pch=15,cex=3)
  points(chrom_dapc$position,chrom_dapc$fst2014int,type='p',col="darkgreen",pch=16,cex=2.5)
  points(chrom_ftempINT$position,chrom_ftempINT$fst2014int,type='p',col="orange",pch=17,cex=2.5)
  #points(chrom_bayes_dapc$position,chrom_bayes_dapc$fst2014int,type='p',col="orange",pch=18,cex=3)

  #add predictor locus for return time that Marine found on Ots12
  #points(20.7,-0.0116, type='p',col="black",pch=18,cex=4)
  
  par(new=T)
  Data.part_2014INT = subset(Data_to_plot_2014INT, subset=chromosome==ijx)
  plot(Data.part_2014INT[,2],Data.part_2014INT[,8],ylim=c(-10,100),xlim=c(0,150),xaxs='i',col="purple",type='l',lwd=2,axes=F,ylab='',xlab='')
  #axis(4,at=c(0,10,20,30,40,50,60,70,80,90,100),labels=FALSE,cex.axis=1.5,cex=1.5)
  #axis(4,at=c(0,20,40,60,80,100),labels=TRUE,cex.axis=2,cex=1.5)
  #legend("topright",legend=c("Temporal Outliers","Bayescan Outliers","DAPC Structure Loci","Moving Average", "# Loci per Window"),pch=c(17,15,16,NA,NA),cex=1.5,pt.cex=2,lwd=2,lty=c(NA,NA,NA,1,1),col=c("green","blue","navy","red","purple"))

############################################################################################################


  data_2002SEG <- data_1998_2002SEG_Fst[data_1998_2002SEG_Fst[,"chromosome"]==ijx,]
  chrom_2002SEG <- just_plot_2002SEG(Nb_bootstrap=1000000, which.chromosome.analysis=ijx ,which.chromosome.plot=ijx ,name="1998_2002SEG")
  #points(chrom_loci$Position,chrom_loci$fst2002seg,type="p",col="black",pch=16,cex=1.2)
  FtempSEG<-ftemp_outliers_SEG[ftemp_outliers_SEG[,"Chrom"]==ijx,]
  points(chrom_bayes$position,chrom_bayes$fst2002seg,type='p',col="blue",pch=15,cex=3)
  points(FtempSEG$position,FtempSEG$fst2002seg,type="p",col="orange",pch=17,cex=2.5)
  points(chrom_dapc$position,chrom_dapc$fst2002seg,type='p',col="darkgreen",pch=16,cex=2.5)
  #points(chrom_bayes_dapc$position,chrom_bayes_dapc$fst2002seg,type='p',col="orange",pch=18,cex=3)
  #bayescan_ftemp_chrom <- bayescan_FtempSEG_overlap[bayescan_FtempSEG_overlap[,"Chrom"]==ijx,]
  #points(bayescan_ftemp_chrom$position,bayescan_ftemp_chrom$fst2002seg,type='p',col="orange",pch=18,cex=4.5)
  
  #add predictor locus for return time that Marine found on Ots12
  #points(20.7,-0.0093, type='p',col="black",pch=18,cex=4)
  
  par(new=T)
  Data.part_2002SEG = subset(Data_to_plot_2002SEG, subset=chromosome==ijx)
  plot(Data.part_2002SEG[,2],Data.part_2002SEG[,8],ylim=c(-10,100),xlim=c(0,150),xaxs='i',col="purple",type='l',lwd=2,axes=F,ylab='',xlab='')
  #axis(4,at=c(0,10,20,30,40,50,60,70,80,90,100),labels=FALSE,cex.axis=1.5,cex=1.5)
  #axis(4,at=c(0,20,40,60,80,100),labels=TRUE,cex.axis=2,cex=1.5)
  
############################################################################################################

  data_2006SEG <- data_1998_2006SEG_Fst[data_1998_2006SEG_Fst[,"chromosome"]==ijx,]
  chrom_2006SEG <- just_plot_2006SEG(Nb_bootstrap=1000000, which.chromosome.analysis=ijx ,which.chromosome.plot=ijx ,name="1998_2006SEG")
  #points(chrom_loci$Position,chrom_loci$fst2006seg,type="p",col="black",pch=16,cex=1.2)
  points(chrom_bayes$position,chrom_bayes$fst2006seg,type='p',col="blue",pch=15,cex=3)
  points(FtempSEG$position,FtempSEG$fst2006seg,type="p",col="orange",pch=17,cex=2.5)
  points(chrom_dapc$position,chrom_dapc$fst2006seg,type='p',col="darkgreen",pch=16,cex=2.5)
  #points(chrom_bayes_dapc$position,chrom_bayes_dapc$fst2006seg,type='p',col="orange",pch=18,cex=3)
  #points(bayescan_ftemp_chrom$position,bayescan_ftemp_chrom$fst2006seg,type='p',col="orange",pch=18,cex=4.5)
  
  #add predictor locus for return time that Marine found on Ots12
  #points(20.7,-0.0035, type='p',col="black",pch=18,cex=4)
  
  par(new=T)
  Data.part_2006SEG = subset(Data_to_plot_2006SEG, subset=chromosome==ijx)
  plot(Data.part_2006SEG[,2],Data.part_2006SEG[,8],ylim=c(-10,100),xlim=c(0,150),xaxs='i',col="purple",type='l',lwd=2,axes=F,ylab='',xlab='')
  #axis(4,at=c(0,10,20,30,40,50,60,70,80,90,100),labels=FALSE,cex.axis=1.5,cex=1.5)
  #axis(4,at=c(0,20,40,60,80,100),labels=TRUE,cex.axis=2,cex=1.5)
  
############################################################################################################

  data_2010SEG <- data_1998_2010SEG_Fst[data_1998_2010SEG_Fst[,"chromosome"]==ijx,]
  chrom_2010SEG <- just_plot_2010SEG(Nb_bootstrap=1000000, which.chromosome.analysis=ijx, which.chromosome.plot=ijx ,name="1998_2010SEG")
  #points(chrom_loci$Position,chrom_loci$fst2010seg,type="p",col="black",pch=16,cex=1.2)
  points(chrom_bayes$position,chrom_bayes$fst2010seg,type='p',col="blue",pch=15,cex=3)
  points(FtempSEG$position,FtempSEG$fst2010seg,type="p",col="orange",pch=17,cex=2.5)
  points(chrom_dapc$position,chrom_dapc$fst2010seg,type='p',col="darkgreen",pch=16,cex=2.5)
  #points(chrom_bayes_dapc$position,chrom_bayes_dapc$fst2010seg,type='p',col="orange",pch=18,cex=3)
  #points(bayescan_ftemp_chrom$position,bayescan_ftemp_chrom$fst2010seg,type='p',col="orange",pch=18,cex=4.5)
    
  #add predictor locus for return time that Marine found on Ots12
  #points(20.7,0.0164, type='p',col="black",pch=18,cex=4)
  
  par(new=T)
  Data.part_2010SEG = subset(Data_to_plot_2010SEG, subset=chromosome==ijx)
  plot(Data.part_2010SEG[,2],Data.part_2010SEG[,8],ylim=c(-10,100),xlim=c(0,150),xaxs='i',col="purple",type='l',lwd=2,axes=F,ylab='',xlab='')
  #axis(4,at=c(0,10,20,30,40,50,60,70,80,90,100),labels=FALSE,cex.axis=1.5,cex=1.5)
  #axis(4,at=c(0,20,40,60,80,100),labels=TRUE,cex.axis=2,cex=1.5)
  

############################################################################################################

  data_2014SEG <- data_1998_2014SEG_Fst[data_1998_2014SEG_Fst[,"chromosome"]==ijx,]
  chrom_2014SEG <- just_plot_2014SEG(Nb_bootstrap=1000000, which.chromosome.analysis=ijx, which.chromosome.plot=ijx ,name="1998_2010SEG")
  #points(chrom_loci$Position,chrom_loci$fst2014seg,type="p",col="black",pch=16,cex=1.2)
  points(chrom_bayes$position,chrom_bayes$fst2014seg,type='p',col="blue",pch=15,cex=3)
  points(FtempSEG$position,FtempSEG$fst2014seg,type="p",col="orange",pch=17,cex=2.5)
  points(chrom_dapc$position,chrom_dapc$fst2014seg,type='p',col="darkgreen",pch=16,cex=2.5)
  #points(chrom_bayes_dapc$position,chrom_bayes_dapc$fst2014seg,type='p',col="orange",pch=18,cex=3)
  #points(bayescan_ftemp_chrom$position,bayescan_ftemp_chrom$fst2014seg,type='p',col="orange",pch=18,cex=4.5)

  #add predictor locus for return time that Marine found on Ots12
  #points(20.7,0.0054, type='p',col="black",pch=18,cex=4)
  
  
  par(new=T)
  Data.part_2014SEG = subset(Data_to_plot_2014SEG, subset=chromosome==ijx)
  plot(Data.part_2014SEG[,2],Data.part_2014SEG[,8],ylim=c(-10,100),xlim=c(0,150),xaxs='i',col="purple",type='l',lwd=2,axes=F,ylab='',xlab='')
  #axis(4,at=c(0,10,20,30,40,50,60,70,80,90,100),labels=FALSE,cex.axis=1.5,cex=1.5)
  #axis(4,at=c(0,20,40,60,80,100),labels=TRUE,cex.axis=2,cex=1.5)



  mtext(paste(ijx),outer=TRUE,line=-1,cex=2,at=0.512)
  mtext("Map Position (cM)",outer=TRUE,side=1,line=2,cex=2,at=0.525)
  mtext("Loci per Window",outer=TRUE, side=4,line=1,cex=2,at=0.525)
  mtext(expression(paste("F"[ST])),outer=TRUE,side=2,line=-2,cex=2,at=0.525)
  

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
#legend("bottom", c("DAPC&Bayes", "Bayes", "DAPC", "Bayes&Temp","Temp", "SW", "Loci per Window"), x.intersp=0.1,text.width=c(0.185,0.17,0.14,0.23,0.24,0.2,0.165),xpd = TRUE, horiz = TRUE, inset = c(0,0), bty = "n", pch = c(18,15,16,18,17,NA,NA), lty=c(NA,NA,NA,NA,NA,1,1),lwd=2, col=c("orange","blue","darkgreen","cyan","green","red","purple"),cex=2.5,pt.cex=c(3,2.5,2.5,3,2.5,2.5,2.5))


#legend for Ots04
legend("bottom", c("Bayes", "Bayes&Temp","Temp", "SW", "Loci per Window"), x.intersp=0.1,text.width=c(0.17,0.23,0.24,0.2,0.165),xpd = TRUE, horiz = TRUE, inset = c(0,0), bty = "n", pch = c(15,18,17,NA,NA), lty=c(NA,NA,NA,1,1),lwd=2, col=c("blue","orange","green","red","purple"),cex=2.5,pt.cex=c(2.5,3,2.5,2.5,2.5))

legend("bottom", c("Bayes", "Temp", "SW"), x.intersp=0.1,text.width=c(0.17,0.24,0.2),xpd = TRUE, horiz = TRUE, inset = c(0,0), bty = "n", pch = c(15,17,NA), lty=c(NA,NA,1),lwd=2, col=c("blue","orange","red"),cex=2.5,pt.cex=c(2.5,2.5,2.5))

  
}

