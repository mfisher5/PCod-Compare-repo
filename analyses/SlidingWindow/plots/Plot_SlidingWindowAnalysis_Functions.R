#################### Plotting Functions for Sliding Window Analysis #################
#
# MF 3/16/2018
#
# This script contains all of the functions required to plot sliding window analyses
# All functions derived from Charlie Waters' code "just_plot_updated_changed_plots"
# To be run AFTER SlidingWindow_MF.R script
# 
#######################################################################################

# Plot just Sliding Window output (best for plotting single chromsomes, full genome) ----------------

just_plot_swa = function(data = swa_output, Nb_bootstrap=100000, Nb_divisions = 150, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="SWA_output_plot")
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
    # set min limit for plot by selecting the min of all data in columns 3:6 of data frame
    min_lim=min(Data.part[,3:6],na.rm=TRUE)-0.02
    # set max limit for plot by selecting the max of all data in columns 3:6 of data frame
    max_lim=max(Data.part[,3:6],na.rm=TRUE)+0.1
    # set max x axis for plot by finding the last marker position
    max_x_lim = max(Data.part[,2])
    # save locations of NAs in 3rd column of data frame ("Fst.Fct")
    where.NA = which(is.na(Data.part[,3]))
    ## plot position of marker on chromosome v. "Fst.Fct" calculated in sla ##
    # if you want to export, create file name
    if( export == TRUE ){
      plotname = paste(name, ijk, sep="_")
      png(paste(plotname, "png", sep="."), width=3000, height=480)
    }
    ## use min/max_lim above to define the y axis
    par(mar=c(5,4,4,2) + 0.1)
    plot(Data.part[,2], Data.part[,3],type="l",ylim=c(min_lim,max_lim), xlim=c(0,max_x_lim), lwd=3, col="red", ylab='',xlab='',cex.lab=1, cex.axis=1)
    # add title, lines to plot
    lines(Data.part[,2], Data.part[,5], col=rgb(0.193,0.205,0.205,0.25))
    lines(Data.part[,2], Data.part[,6], col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)!=0) polygon(c(Data.part[,2][-where.NA], rev(Data.part[,2][-where.NA])), c(Data.part[,5][-where.NA], rev(Data.part[,6][-where.NA])), col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)==0) polygon(c(Data.part[,2], rev(Data.part[,2])), c(Data.part[,5], rev(Data.part[,6])), col=rgb(0.193,0.205,0.205,0.25))    
    abline(h=0)
    mtext(paste(ijk),outer=TRUE,line=-2,cex=2,at=0.515)
    mtext("Map Position (bp)",outer=TRUE,side=1,line=-2,cex=2,at=0.525)
    mtext(expression(paste("F"[ST])),outer=TRUE,side=2,line=-2,cex=2,at=0.525)
    # write out the plot 
    if( export == TRUE ){
      dev.off()
    }
    count = count + 1
    
  }  
  
}



# Plot with Marker Density (best for plotting single chromsomes, full genome) -------------------------

just_plot_md = function(data = swa_output, marker_density = output, Nb_bootstrap=100000, Nb_divisions = 150, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="SWA_output_plot")
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
    # set min limit for plot by selecting the min of all data in columns 3:6 of data frame
    min_lim=min(Data.part[,3:6],na.rm=TRUE)-0.02
    # set max limit for plot by selecting the max of all data in columns 3:6 of data frame
    max_lim=max(Data.part[,3:6],na.rm=TRUE)+0.1
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
    plot(Data.part[,2], Data.part[,3],type="l",ylim=c(min_lim,max_lim), xlim=c(0,max_x_lim), lwd=3, col="red", ylab='',xlab='',cex.lab=1, cex.axis=1)
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
    plot(Data.part[,2],markers_per_window,xlim=c(0,max_x_lim),ylim=c(-10,1.1*max(markers_per_window)),col="black",type='l',lty =2, lwd=1.5, axes = FALSE, ylab='',xlab='')
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








# Plot Multiple Chromosomes with Marker Density (stable axes) ---------------------------------

plot_all = function(data = swa_output, marker_density = output, Nb_bootstrap=100000, Nb_divisions = 150, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="SWA_output_plot")
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






# Plot with Outlier Loci --------------------------------------------------

plot_outliers = function(data = swa_output, outlier_data = outliers, Nb_bootstrap=100000, Nb_divisions = 150, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="SWA_output_outliers_plot")
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
    # subset outlier data frame to only include that chromosome and save as separate data frame
    outliers.part = subset(outlier_data, chromosome %in% ijk)
    print("Number of Outliers on this chromosome:")
    print(length(outliers.part$chromosome))
    # set min limit for plot by selecting the min of all data in columns 3:6 of data frame
    min_lim=min(Data.part[,3:6],na.rm=TRUE)-0.02
    # set max limit for plot by selecting the max of all data in columns 3:6 of data frame
    max_lim=max(max(Data.part[,3:6],na.rm=TRUE),max(outliers.part$fst, na.rm=TRUE))+0.1
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
    plot(Data.part[,2], Data.part[,3],type="l",ylim=c(min_lim,max_lim), xlim=c(0,max_x_lim), lwd=3, col="red", ylab='',xlab='',cex.lab=1, cex.axis=1)
    # add title, lines to plot
    lines(Data.part[,2], Data.part[,5], col=rgb(0.193,0.205,0.205,0.25))
    lines(Data.part[,2], Data.part[,6], col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)!=0) polygon(c(Data.part[,2][-where.NA], rev(Data.part[,2][-where.NA])), c(Data.part[,5][-where.NA], rev(Data.part[,6][-where.NA])), col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)==0) polygon(c(Data.part[,2], rev(Data.part[,2])), c(Data.part[,5], rev(Data.part[,6])), col=rgb(0.193,0.205,0.205,0.25))    
    abline(h=0)
    mtext(paste(ijk),outer=TRUE,line=-2,cex=2,at=0.515)
    mtext("Map Position (bp)",outer=TRUE,side=1,line=-2,cex=2,at=0.525)
    mtext(expression(paste("F"[ST])),outer=TRUE,side=2,line=-2,cex=2,at=0.525)
    # plot over the existing sliding window analysis if there are outliers
    if(length(outliers.part$position) > 1){
      par(new=TRUE)
      plot(outliers.part$position, outliers.part$fst, col="black", pch = 19, ylim=c(min_lim,max_lim), 
           xlim=c(0,max_x_lim),xaxt = "n", yaxt = "n", xlab = "", ylab= "")
    }
    if( export == TRUE ){
      dev.off()
    }
    count = count + 1
    
  }  
  
}










# Plot Two Analyses Stacked --------------------------------------------------
## THIS PLOTTING FUNCTION MUST BE RUN WITHIN SCRIPT `ID_OUTLIER_REGIONS`
just_plot_overlay = function(data1 = swa1_selection, data2 = swa2_selection, Nb_divisions = 150, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="SWA_overlay_plot"){
  if(ncol(data1) < 8){
    print("There is no column 'selection'. Please run ID_Outlier_Regions script.")
  }
  if(which.chromosome.plot != "all"){
    data1 = subset(data1, subset=chromosome %in% which.chromosome.plot)
    data2 = subset(data2, subset=chromosome %in% which.chromosome.plot)
  }
  # set number of chromosomes to the unique values of column chromosome in data (be sure to match column header!!)
  Nb_chromosome = unique(data1$chromosome)
  # set count; used for subsetted the vector "marker density"
  count = 0
  # for each chromosome:
  for (ijk in Nb_chromosome)
  {
    # subset full data frame to only include that chromosome and save as separate data frame
    Data1.part = subset(data1, subset=chromosome %in% ijk)
    Data2.part = subset(data2, subset=chromosome %in% ijk)
    print(head(Data1.part))
    print("Length of moving average data set 1:")
    print(length(Data1.part$position))
    print("Length of moving average data set 2:")
    print(length(Data2.part$position))
    # set max x axis for plot by finding the last marker position
    max_x_lim = max(max(Data1.part$position, na.rm=TRUE), max(Data2.part$position, na.rm=TRUE))
    max_y_lim = max(max(Data1.part$upper_95, na.rm=TRUE), max(Data2.part$upper_95, na.rm=TRUE)) + 0.1
    # save locations of NAs in 3rd column of data frame ("Fst.Fct")
    where.NA.data1 = which(is.na(Data1.part$`Fst/Fct`))
    where.NA.data2 = which(is.na(Data2.part$`Fst/Fct`))
    # save locations potentially under selection in both
    selection1 <- c(Data1.part$position[Data1.part$selection == "Yes"])
    print(length(selection1))
    selection2 <- c(Data2.part$position[Data2.part$selection == "Yes"])
    print(length(selection2))
    ## plot position of marker on chromosome v. "Fst.Fct" calculated in sla ##
    if( export == TRUE ){
      plotname = paste(name, ijk, sep="_")
      png(paste(plotname, "png", sep="."), width=960, height=720)
    }
    ## make the first plot
    par(mfrow=c(2,1), mar=c(0,4,5,2) + 0.1)
    plot(Data1.part$position, Data1.part$`Fst/Fct`,type="l",ylim=c(0,max_y_lim), xlim=c(0,max_x_lim), lwd=3, col="red", ylab='',xaxt = "n", xlab='',cex.lab=1, cex.axis=1)
    # bootstrapping lines to plot
    lines(Data1.part$position, Data1.part$lower_95, col=rgb(0.193,0.205,0.205,0.25))
    lines(Data1.part$position, Data1.part$upper_95, col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA.data1)!=0) polygon(c(Data1.part$position[-where.NA.data1], rev(Data1.part$position[-where.NA.data1])), c(Data1.part$lower_95[-where.NA.data1], rev(Data1.part$upper_95[-where.NA.data1])), col=rgb(0.5,0.205,0.205,0.25))
    if(length(where.NA.data1)==0) polygon(c(Data1.part$position, rev(Data1.part$position)), c(Data1.part$lower_95, rev(Data1.part$upper_95)), col=rgb(0.5,0.205,0.205,0.25))
    abline(h=0)
    abline(v=selection2, col="blue", lty = 3)
    ## make the second plot
    par(mar=c(5,4,0,2) + 0.1)
    plot(Data2.part$position, Data2.part$`Fst/Fct`,type="l",ylim=c(0,max_y_lim), xlim=c(0,max_x_lim), lwd=3, col="blue", ylab='',xlab='',cex.lab=1, cex.axis=1)
    # add title, lines to plot
    lines(Data2.part$position, Data2.part$lower_95, col=rgb(0.193,0.205,0.205,0.25))
    lines(Data2.part$position, Data2.part$upper_95, col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA.data2)!=0) polygon(c(Data2.part$position[-where.NA.data2], rev(Data2.part$position[-where.NA.data2])), c(Data2.part$lower_95[-where.NA.data2], rev(Data2.part$upper_95[-where.NA.data2])), col=rgb(0.193,0.205,0.5,0.25))
    if(length(where.NA.data2)==0) polygon(c(Data2.part$position, rev(Data2.part$position)), c(Data2.part$lower_95, rev(Data2.part$upper_95)), col=rgb(0.193,0.205,0.5,0.25))
    abline(h=0)
    abline(v=selection1, col="red", lty = 3)
    # add axes labels
    mtext("Map Position (bp)",outer=TRUE,side=1,line=-2,cex=2,at=0.525)
    mtext(expression(paste("F"[ST])),outer=TRUE,side=2,line=-2,cex=2,at=0.525)
    mtext(paste(ijk),outer=TRUE,line=-4,cex=2,at=0.515)
    if( export == TRUE ){
      dev.off()
    }
  }
}



# Plot Two Analyses Stacked, Only ID Outlier Regions of Positive Selection --------------------------------------------------
## THIS PLOTTING FUNCTION MUST BE RUN WITHIN SCRIPT `ID_DIVERGENCE_REGIONS`
just_plot_overlay_diverge = function(data1 = swa1_selection, data2 = swa2_selection, Nb_divisions = 150, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="SWA_overlay_plot"){
  if(ncol(data1) < 8){
    print("There is no column 'selection'. Please run ID_Outlier_Regions script.")
  }
  if(which.chromosome.plot != "all"){
    data1 = subset(data1, subset=chromosome %in% which.chromosome.plot)
    data2 = subset(data2, subset=chromosome %in% which.chromosome.plot)
  }
  # set number of chromosomes to the unique values of column chromosome in data (be sure to match column header!!)
  Nb_chromosome = unique(data1$chromosome)
  # set count; used for subsetted the vector "marker density"
  count = 0
  # for each chromosome:
  for (ijk in Nb_chromosome)
  {
    # subset full data frame to only include that chromosome and save as separate data frame
    Data1.part = subset(data1, subset=chromosome %in% ijk)
    Data2.part = subset(data2, subset=chromosome %in% ijk)
    print(head(Data1.part))
    print("Length of moving average data set 1:")
    print(length(Data1.part$position))
    print("Length of moving average data set 2:")
    print(length(Data2.part$position))
    # set max x axis for plot by finding the last marker position
    max_x_lim = max(max(Data1.part$position, na.rm=TRUE), max(Data2.part$position, na.rm=TRUE))
    max_y_lim = max(max(Data1.part$upper_95, na.rm=TRUE), max(Data2.part$upper_95, na.rm=TRUE)) + 0.1
    # save locations of NAs in 3rd column of data frame ("Fst.Fct")
    where.NA.data1 = which(is.na(Data1.part$`Fst/Fct`))
    where.NA.data2 = which(is.na(Data2.part$`Fst/Fct`))
    # save locations potentially under selection in both
    selection1 <- c(Data1.part$position[Data1.part$positive == 1])
    print(length(selection1))
    selection2 <- c(Data2.part$position[Data2.part$positive == 1])
    print(length(selection2))
    ## plot position of marker on chromosome v. "Fst.Fct" calculated in sla ##
    if( export == TRUE ){
      plotname = paste(name, ijk, sep="_")
      png(paste(plotname, "png", sep="."), width=960, height=720)
    }
    ## make the first plot
    par(mfrow=c(2,1), mar=c(0,4,5,2) + 0.1)
    plot(Data1.part$position, Data1.part$`Fst/Fct`,type="l",ylim=c(0,max_y_lim), xlim=c(0,max_x_lim), lwd=3, col="red", ylab='',xaxt = "n", xlab='',cex.lab=1, cex.axis=1)
    # bootstrapping lines to plot
    lines(Data1.part$position, Data1.part$lower_95, col=rgb(0.193,0.205,0.205,0.25))
    lines(Data1.part$position, Data1.part$upper_95, col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA.data1)!=0) polygon(c(Data1.part$position[-where.NA.data1], rev(Data1.part$position[-where.NA.data1])), c(Data1.part$lower_95[-where.NA.data1], rev(Data1.part$upper_95[-where.NA.data1])), col=rgb(0.5,0.205,0.205,0.25))
    if(length(where.NA.data1)==0) polygon(c(Data1.part$position, rev(Data1.part$position)), c(Data1.part$lower_95, rev(Data1.part$upper_95)), col=rgb(0.5,0.205,0.205,0.25))
    abline(h=0)
    abline(v=selection2, col="blue", lty = 3)
    ## make the second plot
    par(mar=c(5,4,0,2) + 0.1)
    plot(Data2.part$position, Data2.part$`Fst/Fct`,type="l",ylim=c(0,max_y_lim), xlim=c(0,max_x_lim), lwd=3, col="blue", ylab='',xlab='',cex.lab=1, cex.axis=1)
    # add title, lines to plot
    lines(Data2.part$position, Data2.part$lower_95, col=rgb(0.193,0.205,0.205,0.25))
    lines(Data2.part$position, Data2.part$upper_95, col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA.data2)!=0) polygon(c(Data2.part$position[-where.NA.data2], rev(Data2.part$position[-where.NA.data2])), c(Data2.part$lower_95[-where.NA.data2], rev(Data2.part$upper_95[-where.NA.data2])), col=rgb(0.193,0.205,0.5,0.25))
    if(length(where.NA.data2)==0) polygon(c(Data2.part$position, rev(Data2.part$position)), c(Data2.part$lower_95, rev(Data2.part$upper_95)), col=rgb(0.193,0.205,0.5,0.25))
    abline(h=0)
    abline(v=selection1, col="red", lty = 3)
    # add axes labels
    mtext("Map Position (bp)",outer=TRUE,side=1,line=-2,cex=2,at=0.525)
    mtext(expression(paste("F"[ST])),outer=TRUE,side=2,line=-2,cex=2,at=0.525)
    mtext(paste(ijk),outer=TRUE,line=-4,cex=2,at=0.515)
    if( export == TRUE ){
      dev.off()
    }
  }
}



# Plot Sliding Window Analyses, Only ID Outlier Regions of Positive Selection --------------------------------------------------
## THIS PLOTTING FUNCTION MUST BE RUN WITHIN SCRIPT `ID_DIVERGENCE_REGIONS`
just_plot_diverge = function(data1 = swa_selection, Nb_divisions = 150, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="SWA_divergence_plot"){
  if(ncol(data1) < 8){
    print("There is no column 'selection'. Please run ID_Outlier_Regions script.")
  }
  if(which.chromosome.plot != "all"){
    data1 = subset(data1, subset=chromosome %in% which.chromosome.plot)
  }
  # set number of chromosomes to the unique values of column chromosome in data (be sure to match column header!!)
  Nb_chromosome = unique(data1$chromosome)
  # set count; used for subsetted the vector "marker density"
  count = 0
  # for each chromosome:
  for (ijk in Nb_chromosome)
  {
    # subset full data frame to only include that chromosome and save as separate data frame
    Data1.part = subset(data1, subset=chromosome %in% ijk)
    print(head(Data1.part))
    print("Length of moving average data set 1:")
    print(length(Data1.part$position))
    # set max x axis for plot by finding the last marker position
    max_x_lim = max(max(Data1.part$position, na.rm=TRUE))
    max_y_lim = max(max(Data1.part$upper_95, na.rm=TRUE)) + 0.1
    # save locations of NAs in 3rd column of data frame ("Fst.Fct")
    where.NA.data1 = which(is.na(Data1.part$`Fst/Fct`))
    # save locations potentially under selection in both
    selection1 <- c(Data1.part$position[Data1.part$positive == 1])
    print(length(selection1))
    ## plot position of marker on chromosome v. "Fst.Fct" calculated in sla ##
    if( export == TRUE ){
      plotname = paste(name, ijk, sep="_")
      png(paste(plotname, "png", sep="."), width=960, height=480)
    }
    ## make the first plot
    par(mar=c(5,4,4,2) + 0.1)
    plot(Data1.part$position, Data1.part$`Fst/Fct`,type="l",ylim=c(0,max_y_lim), xlim=c(0,max_x_lim), lwd=3, col="red", ylab='',xaxt = "n", xlab='',cex.lab=1, cex.axis=1)
    # bootstrapping lines to plot
    lines(Data1.part$position, Data1.part$lower_95, col=rgb(0.193,0.205,0.205,0.25))
    lines(Data1.part$position, Data1.part$upper_95, col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA.data1)!=0) polygon(c(Data1.part$position[-where.NA.data1], rev(Data1.part$position[-where.NA.data1])), c(Data1.part$lower_95[-where.NA.data1], rev(Data1.part$upper_95[-where.NA.data1])), col=rgb(0.5,0.205,0.205,0.25))
    if(length(where.NA.data1)==0) polygon(c(Data1.part$position, rev(Data1.part$position)), c(Data1.part$lower_95, rev(Data1.part$upper_95)), col=rgb(0.5,0.205,0.205,0.25))
    abline(h=0)
    abline(v=selection1, col="red", lty = 3)
    # add axes labels
    mtext("Map Position (bp)",outer=TRUE,side=1,line=-2,cex=2,at=0.525)
    mtext(expression(paste("F"[ST])),outer=TRUE,side=2,line=-2,cex=2,at=0.525)
    mtext(paste(ijk),outer=TRUE,line=-4,cex=2,at=0.515)
    if( export == TRUE ){
      dev.off()
    }
  }
}



# Plot Heterozygosity Sliding Window --------------------------------------------------
## THIS PLOTTING FUNCTION MUST BE RUN WITHIN SCRIPT `SLIDINGWINDOW_HO_EAST/WEST_MF.R`

plot_het = function(data = swa_het_output, Nb_bootstrap=100000, Nb_divisions = 150, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="SWA_Het_output")
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
    # set max x axis for plot by finding the last marker position
    max_x_lim = max(Data.part[,2])
    # save locations of NAs in 3rd column of data frame ("Het")
    where.NA = which(is.na(Data.part[,3]))
    ## plot position of marker on chromosome v. "Het" calculated in sla ##
    # if you want to export, create file name
    if( export == TRUE ){
      plotname = paste(name, ijk, sep="_")
      png(paste(plotname, "png", sep="."), width=960, height=480)
    }
    ## use min/max_lim above to define the y axis
    par(mar=c(5,4,4,2) + 0.1)
    plot(Data.part[,2], Data.part[,3],type="l",ylim=c(0,1.0), xlim=c(0,max_x_lim), lwd=3, col="red", ylab='',xlab='',cex.lab=1, cex.axis=1)
    # add title, lines to plot
    abline(h=0)
    mtext(paste(ijk),outer=TRUE,line=-2,cex=2,at=0.515)
    mtext("Map Position (bp)",outer=TRUE,side=1,line=-2,cex=2,at=0.525)
    mtext(expression(paste("H"[observed])),outer=TRUE,side=2,line=-2,cex=2,at=0.525)
    # write out the plot 
    if( export == TRUE ){
      dev.off()
    }
    count = count + 1
    
  }  
  
}

# Plot Two Heterozygosity Sliding Windows --------------------------------------------------
## THIS PLOTTING FUNCTION MUST BE RUN WITHIN SCRIPT `PLOT_HETEROZYGOSITY_SLIDINGWINDOW.R`

plot_het_overlay = function(data1 = swa_het_output_east, data2 = swa_het_output_west, legend.text.data1 = "East", legend.text.data2 = "West", Nb_bootstrap=100000, Nb_divisions = 150, which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="SWA_Het_East_West_overlay")
{
  # subset the data if not looking at all chromosomes; otherwise, load in all data
  if(which.chromosome.analysis!="all") data1 = subset(data1, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis!="all") data2 = subset(data2, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis=="all") data1 = data1
  if(which.chromosome.analysis=="all") data2 = data2
  # set number of chromosomes to the unique values of column chromosome in data (be sure to match column header!!)
  Nb_chromosome = unique(data1$chromosome)
  # set count; used for subsetted the vector "marker density"
  count = 0
  # for each chromosome:
  for (ijk in Nb_chromosome)
  {
    # subset full data frame to only include that chromosome and save as separate data frame
    Data1.part = subset(data1, subset=chromosome %in% ijk)
    Data2.part = subset(data2, subset=chromosome %in% ijk)
    print(head(Data1.part))
    print("Length of moving average data 1:")
    print(length(Data1.part[,2]))
    print(head(Data2.part))
    print("Length of moving average data 2:")
    print(length(Data2.part[,2]))
    # set max x axis for plot by finding the last marker position
    max_x_lim = max(max(Data1.part[,2]), max(Data2.part[,2]))
    # save locations of NAs in 3rd column of data frame ("Het")
    where.NA1 = which(is.na(Data1.part[,3]))
    where.NA2 = which(is.na(Data2.part[,3]))
    ## plot position of marker on chromosome v. "Het" calculated in sla ##
    # if you want to export, create file name
    if( export == TRUE ){
      plotname = paste(name, ijk, sep="_")
      png(paste(plotname, "png", sep="."), width=960, height=480)
    }
    ## plot first data set. use min/max_lim above to define the x axis
    par(mar=c(5,4,4,2) + 0.1)
    plot(Data1.part[,2], Data1.part[,3],type="l",ylim=c(0,0.5), xlim=c(0,max_x_lim), lwd=3, col="red", ylab='',xlab='',cex.lab=1, cex.axis=1)
    ## plot second data set. use min/max_lim above to define the x axis
    par(new=TRUE)
    plot(Data2.part[,2], Data2.part[,3],type="l",ylim=c(0,0.5), xlim=c(0,max_x_lim), lwd=3, col="blue", ylab='',xlab='',xaxt = "n", yaxt = "n")
    # add title, lines to plot
    abline(h=0)
    mtext(paste(ijk),outer=TRUE,line=-2,cex=2,at=0.515)
    mtext("Map Position (bp)",outer=TRUE,side=1,line=-2,cex=2,at=0.525)
    mtext(expression(paste("H"[observed])),outer=TRUE,side=2,line=-2,cex=2,at=0.525)
    legend("topright", legend=c(legend.text.data1, legend.text.data2), fill = c("red", "blue"))
    # write out the plot 
    if( export == TRUE ){
      dev.off()
    }
    count = count + 1
    
  }  
  
}

# Plot Three Analyses, Lines only ---------------------------------

plot_lines_overlay = function(data1 = swa_output, data2 = swa_output2, data3 = swa_output3, Nb_divisions = 150, legend.text= c("Analysis 1","Analysis 2", "Analysis 3"), which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="SWA_lines_overlay")
{
  # subset the data if not looking at all chromosomes; otherwise, load in all data
  if(which.chromosome.analysis!="all") data1 = subset(data1, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis!="all") data2 = subset(data2, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis!="all") data3 = subset(data3, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis=="all") data1 = data1
  if(which.chromosome.analysis=="all") data2 = data2
  if(which.chromosome.analysis=="all") data3 = data3
  # set number of chromosomes to the unique values of column chromosome in data (be sure to match column header!!)
  Nb_chromosome = unique(data1$chromosome)
  # set count; used for subsetted the vector "marker density"
  count = 0
  # for each chromosome:
  for (ijk in Nb_chromosome)
  {
    # subset full data frame to only include that chromosome and save as separate data frame
    Data1.part = subset(data1, subset=chromosome %in% ijk)
    Data2.part = subset(data2, subset=chromosome %in% ijk)
    Data3.part = subset(data3, subset=chromosome %in% ijk)
    print(head(Data1.part))
    print("Length of moving average data 1:")
    print(length(Data1.part[,2]))
    print(head(Data2.part))
    print("Length of moving average data 2:")
    print(length(Data2.part[,2]))
    print(head(Data3.part))
    print("Length of moving average data 3:")
    print(length(Data3.part[,2]))
    # set max x axis for plot by finding the last marker position
    max_x_lim = max(max(Data1.part[,2]), max(Data2.part[,2]), max(Data3.part[,2]))
    # save locations of NAs in 3rd column of data frame ("Het")
    where.NA1 = which(is.na(Data1.part[,3]))
    where.NA2 = which(is.na(Data2.part[,3]))
    where.NA3 = which(is.na(Data3.part[,3]))
    
    ## plot position of marker on chromosome v. "Het" calculated in sla ##
    # if you want to export, create file name
    if( export == TRUE ){
      plotname = paste(name, ijk, sep="_")
      png(paste(plotname, "png", sep="."), width=960, height=480)
    }
    ## plot first data set. use min/max_lim above to define the x axis
    par(mar=c(5,4,4,2) + 0.1)
    plot(Data1.part[,2], Data1.part[,3],type="l",ylim=c(0,1), xlim=c(0,max_x_lim), lwd=3, col="black", ylab='',xlab='',cex.lab=1, cex.axis=1)
    ## plot second data set. use min/max_lim above to define the x axis
    par(new=TRUE)
    plot(Data2.part[,2], Data2.part[,3],type="l",ylim=c(0,1), xlim=c(0,max_x_lim), lwd=3, col="red", ylab='',xlab='',xaxt = "n", yaxt = "n")
    ## plot third data set. use min/max_lim above to define the x axis
    par(new=TRUE)
    plot(Data3.part[,2], Data3.part[,3],type="l",ylim=c(0,1), xlim=c(0,max_x_lim), lwd=3, col="blue", ylab='',xlab='',xaxt = "n", yaxt = "n")
    # add title, lines to plot
    abline(h=0)
    mtext(paste(ijk),outer=TRUE,line=-2,cex=2,at=0.515)
    mtext("Map Position (bp)",outer=TRUE,side=1,line=-2,cex=2,at=0.525)
    mtext(expression(paste("F"[st])),outer=TRUE,side=2,line=-2,cex=2,at=0.525)
    legend("topright", legend=c(legend.text), fill = c("black", "red", "blue"))
    # write out the plot 
    if( export == TRUE ){
      dev.off()
    }
    count = count + 1
    
  }  
  
}





# Plot Three Analyses, Lines only, on one plot w/ Islands ---------------------------------

plot_lines_overlay_divergence = function(data1 = swa_output_selection, data2 = swa_output2_selection, data3 = swa_output3_selection, Nb_divisions = 150, legend.text= c("Analysis 1","Analysis 2", "Analysis 3"), which.chromosome.analysis="all", which.chromosome.plot="all",export = TRUE, name="SWA_lines_overlay")
{
  # subset the data if not looking at all chromosomes; otherwise, load in all data
  if(which.chromosome.analysis!="all") data1 = subset(data1, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis!="all") data2 = subset(data2, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis!="all") data3 = subset(data3, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis=="all") data1 = data1
  if(which.chromosome.analysis=="all") data2 = data2
  if(which.chromosome.analysis=="all") data3 = data3
  # set number of chromosomes to the unique values of column chromosome in data (be sure to match column header!!)
  Nb_chromosome = unique(data1$chromosome)
  # set count; used for subsetted the vector "marker density"
  count = 0
  # for each chromosome:
  for (ijk in Nb_chromosome)
  {
    # subset full data frame to only include that chromosome and save as separate data frame
    Data1.part = subset(data1, subset=chromosome %in% ijk)
    Data2.part = subset(data2, subset=chromosome %in% ijk)
    Data3.part = subset(data3, subset=chromosome %in% ijk)
    print(head(Data1.part))
    print("Length of moving average data 1:")
    print(length(Data1.part[,2]))
    print(head(Data2.part))
    print("Length of moving average data 2:")
    print(length(Data2.part[,2]))
    print(head(Data3.part))
    print("Length of moving average data 3:")
    print(length(Data3.part[,2]))
    # set max x axis for plot by finding the last marker position
    max_x_lim = max(max(Data1.part[,2]), max(Data2.part[,2]), max(Data3.part[,2]))
    # save locations of NAs in 3rd column of data frame ("Het")
    where.NA1 = which(is.na(Data1.part[,3]))
    where.NA2 = which(is.na(Data2.part[,3]))
    where.NA3 = which(is.na(Data3.part[,3]))
    
    # save locations potentially under selection in both
    print("Number of Windows Under Selection in Data sets:")
    selection1 <- c(Data1.part$position[Data1.part$positive == 1])
    print(length(selection1))
    selection2 <- c(Data2.part$position[Data2.part$positive == 1])
    print(length(selection2))
    selection3 <- c(Data3.part$position[Data3.part$positive == 1])
    print(length(selection3))
    
    ## plot position of marker on chromosome v. "Het" calculated in sla ##
    # if you want to export, create file name
    if( export == TRUE ){
      plotname = paste(name, ijk, sep="_")
      png(paste(plotname, "png", sep="."), width=960, height=480)
    }
    ## plot first data set. use min/max_lim above to define the x axis
    par(mar=c(5,4,4,2) + 0.1)
    plot(Data1.part[,2], Data1.part[,3],type="l",ylim=c(0,1), xlim=c(0,max_x_lim), lwd=3, col="black", ylab='',xlab='',cex.lab=1, cex.axis=1)
    abline(v=selection1, col="black", lty = 3)
    ## plot second data set. use min/max_lim above to define the x axis
    par(new=TRUE)
    plot(Data2.part[,2], Data2.part[,3],type="l",ylim=c(0,1), xlim=c(0,max_x_lim), lwd=3, col="red", ylab='',xlab='',xaxt = "n", yaxt = "n")
    abline(v=selection2, col="red", lty = 2)
    ## plot third data set. use min/max_lim above to define the x axis
    par(new=TRUE)
    plot(Data3.part[,2], Data3.part[,3],type="l",ylim=c(0,1), xlim=c(0,max_x_lim), lwd=3, col="blue", ylab='',xlab='',xaxt = "n", yaxt = "n")
    abline(v=selection3, col="blue", lty = 4)
    # add title, lines to plot
    abline(h=0)
    mtext(paste(ijk),outer=TRUE,line=-2,cex=2,at=0.515)
    mtext("Map Position (bp)",outer=TRUE,side=1,line=-2,cex=2,at=0.525)
    mtext(expression(paste("F"[st])),outer=TRUE,side=2,line=-2,cex=2,at=0.525)
    legend("topright", legend=c(legend.text), fill = c("black", "red", "blue"))
    # write out the plot 
    if( export == TRUE ){
      dev.off()
    }
    count = count + 1
    
  }  
  
}








