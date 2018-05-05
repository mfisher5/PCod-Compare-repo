################### Plot Sliding Window Analysis with All Loci, Highlighting Outliers #####################
#
# This script will overlay per locus FST on sliding window analysis output
#     and highlight which loci are outliers. 
# The code for the plotting function is included in this script.
#
# Mary Fisher 5/3/2018, based on plotting script from Charlie Waters
#
###########################################################################################################


# Load Packages -----------------------------------------------------------
library(dplyr)



# Read in Data ------------------------------------------------------------
setwd("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/SlidingWindow")

## sliding window analysis output file. Use filtered version
swa_output = read.table("West/batch_8_final_filtered_west_2reg_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150_FILTERED.txt", header = TRUE, sep = "\t")
dim(swa_output)
head(swa_output) ## heading should be:  chromosome,position,Fst.Fct,Mean_boostrap,lower_95,upper_95,pvalue


## sliding window analysis input file (will provide FST per locus)
swa_input <- read.table("West/batch_8_SWA_input_west_sorted.txt", header = TRUE, sep="\t", colClasses = c("character", "numeric", "character", "numeric"))
head(swa_input) ## heading should be:  Locus,fst,chromosome,position


## list of outlier loci to highlight. should have locus names in first column; adjust "colClasses" per your input file
outlierfile = read.table("../outliers/batch_8_final_filtered_aligned_WEST_outliers_2progs_forSWA.txt", header = TRUE, colClasses = c("character"))
head(outlierfile)



# Prepare Data for Plotting Function --------------------------------------
## get list of outlier loci from data frame
outlier_loci <- outlierfile$Locus
str(outlier_loci) # should be characters



# Plotting Function -------------------------------------------------------
## Key for editing plotting function
##    1. Change output picture size: 98
##    2. Change plotting parameters for average line of sliding window analysis: 102
##    3. Change plotting parameters for CI polygons: 106-107
##    4. Change plotting parameters for CI upper / lower limit lines: 104-105
##    5. Change plotting parameters for loci: 114
##    6. Change plotting parameters for outlier loci: 119

plot_loci_overlay = function(data = swa_output, loci_data = swa_input, outlier_data = outliers, Nb_bootstrap=100000, Nb_divisions = 150, color="red", ylim = NULL, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="SWA_output_outliers_plot")
{
  # subset the data if not looking at all chromosomes; otherwise, load in all data
  if(which.chromosome.analysis!="all") data = subset(data, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis=="all") data = data
  # set number of chromosomes to the unique values of column chromosome in data (be sure to match column header!!)
  Nb_chromosome = unique(data$chromosome)
  # set count; used for subsetted the vector "marker density"
  count = 0
  # for each chromosome:
  for (ijk in Nb_chromosome)
  {
    # subset full data frame to only include that chromosome and save as separate data frame
    Data.part = subset(data, subset=chromosome %in% ijk)
    print(head(Data.part))
    print("Length of moving average data:")
    print(length(Data.part[,2]))
    # subset loci data frame to only include that chromosome and save as separate data frame
    loci.part = subset(loci_data, chromosome %in% ijk)
    print("Number of Loci on this chromosome:")
    print(length(loci.part$chromosome))
    # subset loci data frame to only include outliers
    outliers.part = subset(loci.part, Locus %in% outlier_data)
    print("Number of Outliers on this chromosome:")
    print(length(outliers.part$chromosome))
    # set min limit for plot by selecting the min of all data in columns 3:6 of data frame
    min_lim=min(Data.part[,3:6],na.rm=TRUE)-0.02
    # set max limit for plot by selecting the max of all data in columns 3:6 of data frame
    if(ylim == "NULL"){
      max_lim=max(max(Data.part[,3:6],na.rm=TRUE),max(loci.part$fst, na.rm=TRUE))+0.1
    } else {max_lim = ylim}
    # set max x axis for plot by finding the last marker position
    max_x_lim = max(Data.part[,2])
    # save locations of NAs in 3rd column of data frame ("Fst.Fct")
    where.NA = which(is.na(Data.part[,3]))
    # subset marker data based on chromosome you are looking at
    #marker_density_subset = subset(marker_density, chrom %in% ijk)
    #markers_per_window = marker_density_subset$num_markers
    #print(head(marker_density_subset))
    #print("length of marker density vector:")
    #print(length(markers_per_window))
    ## plot position of marker on chromosome v. "Fst.Fct" calculated in sla ##
    # if you want to export, create file name
    if( export == TRUE ){
      plotname = paste(name, ijk, sep="_")
      png(paste(plotname, "png", sep="."), width=960, height=480)
    }
    ## use min/max_lim above to define the y axis
    par(mar=c(5,4,4,2) + 0.1)
    plot(Data.part[,2], Data.part[,3],type="l",ylim=c(min_lim,max_lim), xlim=c(0,max_x_lim), lwd=3, col=color, ylab='',xlab='',cex.lab=1, cex.axis=1)
    # add title, lines to plot
    if(length(where.NA)!=0) lines(Data.part[,2][-where.NA], Data.part[,5][-where.NA], col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)!=0) lines(Data.part[,2][-where.NA], Data.part[,6][-where.NA], col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)!=0) polygon(c(Data.part[,2][-where.NA], rev(Data.part[,2][-where.NA])), c(Data.part[,5][-where.NA], rev(Data.part[,6][-where.NA])), col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)==0) polygon(c(Data.part[,2], rev(Data.part[,2])), c(Data.part[,5], rev(Data.part[,6])), col=rgb(0.193,0.205,0.205,0.25))    
    abline(h=0)
    mtext(paste(ijk),outer=TRUE,line=-2,cex=2,at=0.515)
    mtext("Map Position (bp)",outer=TRUE,side=1,line=-2,cex=2,at=0.525)
    mtext(expression(paste("F"[ST])),outer=TRUE,side=2,line=-2,cex=2,at=0.525)
    # plot loci over the existing sliding window analysis
    par(new=TRUE)
    plot(loci.part$position, loci.part$fst, col = "black", pch = 15, ylim = c(min_lim,max_lim),
         xlim = c(0, max_x_lim), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    # change color of outlier markers by plotting over loci markers
    if(length(outliers.part$position) > 0){
      par(new=TRUE)
      print(outliers.part)
      plot(outliers.part$position, outliers.part$fst, col=color, pch = 15, ylim=c(min_lim,max_lim), 
           xlim=c(0,max_x_lim),xaxt = "n", yaxt = "n", xlab = "", ylab= "")
    }
    if( export == TRUE ){
      dev.off()
    }
    count = count + 1
    
  }  
  print(paste("Finished plotting ",count,"chromosomes"))
  
}













# Plot Data ---------------------------------------------------------------

plot_loci_overlay(data = swa_output, loci_data = swa_input, outlier_data = outlier_loci, Nb_bootstrap=100000, Nb_divisions = 150, color="red", which.chromosome.analysis="all", which.chromosome.plot="all",ylim = 1.0, export = TRUE, name="plots/EastWest/EastWest_SWA_AllLoci_plot")

plot_loci_overlay(data = swa_output, loci_data = swa_input, outlier_data = outlier_loci, Nb_bootstrap=100000, Nb_divisions = 150, color="deepskyblue4", which.chromosome.analysis="all", which.chromosome.plot="all",ylim = 0.5, export = TRUE, name="plots/West/West_SWA_AllLoci_plot")

plot_loci_overlay(data = swa_output, loci_data = swa_input, outlier_data = outlier_loci, Nb_bootstrap=100000, Nb_divisions = 150, color="darkorchid2", which.chromosome.analysis="all", which.chromosome.plot="all",ylim = 0.25, export = FALSE, name="plots/East/East_SWA_AllLoci_plot")


