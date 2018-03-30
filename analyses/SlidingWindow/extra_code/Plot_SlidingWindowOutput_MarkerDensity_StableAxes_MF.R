############# Plotting Sliding Window Analysis / Marker Density, Stable Axis #########
#
# Adopted from Charlie Waters' script "just_plot_updated_changed_plots"
# To be run AFTER SlidingWindow_MF.R script
# Allows for the use of plotting function across ALL chromosomes as input
#
# MF edited for Pacific cod 2/26/2018
#
#########################################################################


# Install Packages --------------------------------------------------------

install.packages("ggplot2")
library(ggplot2)
library(dplyr)


# Load data ---------------------------------------------------------------
## set working directory
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/SlidingWindow")

## input data is the same as for sliding window analysis (the sorted alignment table output as a file in previous script)
marker_data = read.table("batch_8_final_filtered_SWA_input_eastwest_sorted.txt", header = TRUE, sep = "\t")
View(marker_data)



# Set Parameters ----------------------------------------------------------
## these should match the parameters used in sliding window analysis
window_size=250000
divisions=150
Num_marker_data = nrow(marker_data)
Num_chromosome = unique(marker_data$chromosome)
chromosomes_all= rep(Num_chromosome,each=150) # `each` argument should == divisions
Positions_all=c()
markers_per_window=c()
par(mfrow=c(1,1), mar=c(2,2,2,2))



# Calculate # Markers per Window ------------------------------------------
## CHARLIE'S CODE: for loop to calculate markers per window on each chromosome
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


## save output as a matrix
output=c()   ##be sure to clear output from previous runs
output=data.frame(chromosomes_all,Positions_all,markers_per_window)
colnames(output)=c("chrom","position","num_markers")
View(output)


## check calculated number markers per window for first window of chromosomes 1 - 3
##-- input names of first three chromosomes here
chr1 <- "LG01"
chr2 <- "LG02"
chr3 <- "LG03"
sigma <- window_size

##-- set full window size
full_window = sigma*6; full_window

##-- first chromosome
markers_in_window <- markers_per_window[1]
marker_data_chr1 <- marker_data[marker_data$chromosome == chr1,]
end_window <- marker_data_chr1$position[1] + full_window
counted_markers_in_window <- nrow(marker_data_chr1[marker_data_chr1$position < end_window,])
if(markers_in_window == counted_markers_in_window){
  print("FUNCTION VERIFIED")
} else{
  print("ERROR:")
  print(markers_in_window)
  print(counted_markers_in_window)
}

##-- second chromosome
markers_in_window <- markers_per_window[151]
marker_data_chr2 <- marker_data[marker_data$chromosome == chr2,]
end_window <- marker_data_chr2$position[1] + full_window
counted_markers_in_window <- nrow(marker_data_chr2[marker_data_chr2$position < end_window,])
if(markers_in_window == counted_markers_in_window){
  print("FUNCTION VERIFIED")
} else{
  print("ERROR")
  print(markers_in_window)
  print(counted_markers_in_window)
}

##-- third chromosome
markers_in_window <- markers_per_window[301]
marker_data_chr3 <- marker_data[marker_data$chromosome == chr3,]
end_window <- marker_data_chr3$position[1] + full_window
counted_markers_in_window <- nrow(marker_data_chr3[marker_data_chr3$position < end_window,])
if(markers_in_window == counted_markers_in_window){
  print("FUNCTION VERIFIED")
} else{
  print("ERROR")
  print(markers_in_window)
  print(counted_markers_in_window)
}




# Explore Number Markers per Window ------------------------------------------

## find average number of markers per window
mean(markers_per_window)

## plot number markers per window
qplot(markers_per_window, geom="histogram",
      binwidth = 1,
      main = "Markers per Sliding Window, East v. West\nAcross All Linkage Groups",
      xlab = "# Markers", 
      ylab = "# Windows") +
  annotate("text", x = 25, y = 300, label="Mean:\n10.97")


# Plot Markers per Window over SWA ----------------------------------------
## Load in output file from Kot's function to calculate sliding window
sla_output = read.table("batch_8_final_filtered_eastwest_output_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150.txt", header = TRUE, sep = "\t")
View(sla_output)
dim(sla_output)


## READ IN CHARLIE'S PLOTTING FUNCTION (adjusted)
just_plot = function(data = swa_output, marker_density = output, Nb_bootstrap=100000, Nb_divisions = 150, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="SWA_output_plot")
{
  # subset the data if not looking at all chromosomes; otherwise, load in all data
  if(which.chromosome.analysis!="all") data = subset(data, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis=="all") data = data
  # set number of chromosomes to the unique values of column chromosome in data (be sure to match column header!!)
  Nb_chromosome = unique(data$chromosome)
  # set count; used for subsetted the vector "marker density"
  count = 0
  # select axis limit for marker density; will apply to all plots
  markers_ylim <- 1.1*max(marker_density$num_markers)
  # for each chromosome:
  for (ijk in Nb_chromosome)
  {
    # subset full data frame to only include that chromosome and save as separate data frame
    Data.part = subset(data, subset=chromosome %in% ijk)
    print(head(Data.part))
    print("Length of moving average data:")
    print(length(Data.part[,2]))
    # set max x axis for plot by finding the last marker position
    max_x_lim = max(Data.part[,2])
    # save locations of NAs in 3rd column of data frame ("Fst.Fct")
    where.NA = which(is.na(Data.part[,3]))
    # subset marker data based on chromosome you are looking at
    marker_density_subset = subset(marker_density, chrom %in% ijk)
    markers_per_window = marker_density_subset$num_markers
    print(head(marker_density_subset))
    print("length of marker density vector:")
    print(length(markers_per_window))
    ## plot position of marker on chromosome v. "Fst.Fct" calculated in sla ##
    # if you want to export, create file name
    if( export == TRUE ){
      plotname = paste(name, ijk, sep="_")
      png(paste(plotname, "png", sep="."), width=960, height=480)
    }
    ## use min/max_lim above to define the y axis
    par(mar=c(5,4,4,2) + 0.1)
    plot(Data.part[,2], Data.part[,3],type="l",ylim=c(0,1.0), xlim=c(0,max_x_lim), lwd=3, col="red", ylab='',xlab='',cex.lab=1, cex.axis=1)
    # add title, lines to plot
    lines(Data.part[,2], Data.part[,5], col=rgb(0.193,0.205,0.205,0.25))
    lines(Data.part[,2], Data.part[,6], col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)!=0) polygon(c(Data.part[,2][-where.NA], rev(Data.part[,2][-where.NA])), c(Data.part[,5][-where.NA], rev(Data.part[,6][-where.NA])), col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)==0) polygon(c(Data.part[,2], rev(Data.part[,2])), c(Data.part[,5], rev(Data.part[,6])), col=rgb(0.193,0.205,0.205,0.25))    
    abline(h=0)
    mtext(paste(ijk),outer=TRUE,line=-2,cex=2,at=0.515)
    mtext("Map Position (bp)",outer=TRUE,side=1,line=-2,cex=2,at=0.525)
    mtext("Loci per Window",outer=TRUE, side=4,line=-2,cex=1,at=0.25)
    mtext(expression(paste("F"[ST])),outer=TRUE,side=2,line=-2,cex=2,at=0.525)
    # plot over the existing sliding window analysis:
    par(new=TRUE)
    # plot # markers per window (x axis stays the same) in dashed line (y=0 starts a bit higher for easier visualization)
    plot(Data.part[,2],markers_per_window,xlim=c(0,max_x_lim),ylim=c(-10,markers_ylim),col="black",type='l',lty =2, lwd=1.5, axes = FALSE, ylab='',xlab='')
    abline(h=0, lty=2, col="black")
    # add tic marks on second y axis, with labels
    axis(4,at=c(0,5,10,15,20,25,30,35,40,45,50),labels=TRUE,cex.axis=1,cex=1)
    # add y axis label
    mtext("Loci per Window",outer = TRUE, side=4,line=1,cex=1.5, at=0.525)
    # write out the plot 
    if( export == TRUE ){
      dev.off()
    }
    count = count + 1
    
  }  
  
}



## Run plotting function on all chromosomes
##-- double check dimensions
if(length(sla_output$position) != length(markers_per_window)){
  print("ERROR: position and marker density vectors are not of equal lengths. you cannot run the just-plot function until you subset markers_per_window.")
} else{print("continue to plotting function.")}

##-- plot
just_plot(data = sla_output, marker_density = output, Nb_bootstrap=100000, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="batch_8_eastwest_SWA")
?par





## Run plotting function on specific chromosomes

##-- subset according to # divisions per chromosome
markers_per_window_lg01 <- markers_per_window[1:150]
##-- check dimensions
if(length(sla_output$position) != length(markers_per_window_lg01)){
  print("ERROR: position and marker density vectors are not of equal lengths. you cannot run the just-plot function until you subset markers_per_window.")
} else{print("continue to plotting function.")}
##-- plot
just_plot(data = sla_output, marker_density = output, Nb_bootstrap=100000, which.chromosome.analysis="all", which.chromosome.plot="all",name="batch_8_eastwest_output")
?par


