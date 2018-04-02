plot_fst_het = function(data = swa_output,
                        het1 = swa_het_output_east,
                        het2 = swa_het_output_west,
                        legend.text.het1 = "East", legend.text.het2 = "West",
                        Nb_bootstrap_fst=100000, Nb_divisions_fst = 150, 
                        Nb_bootstrap_het=100000, Nb_divisions_het = 150,
                        which.chromosome.analysis="all", which.chromosome.plot="all",
                        export = TRUE, name="SWA_output_plot") {
  # subset the data if not looking at all chromosomes; otherwise, load in all data
  if(which.chromosome.analysis!="all") het1 = subset(het1, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis!="all") het2 = subset(het2, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis!="all") data = subset(data, subset=chromosome %in% which.chromosome.analysis)
  if(which.chromosome.analysis=="all") data = data
  if(which.chromosome.analysis=="all") het1 = het1
  if(which.chromosome.analysis=="all") het2 = het2
  # set number of chromosomes to the unique values of column chromosome in data (be sure to match column header!!)
  Nb_chromosome = unique(data$chromosome)
  # set count; used for subsetted the vector "marker density"
  count = 0
  # select y axis limit for heterozygosity; will apply to all plots
  het_ylim <- 1.1*max(het1[,3])
  # select y axis limit for fst; will apply to all plots
  fst_ylim <- 1.1*max(data[,5])
  # for each chromosome:
  for (ijk in Nb_chromosome)
  {
    # subset fst data to only include that chromosome and save as separate data frame
    Data.part = subset(data, subset=chromosome %in% ijk)
    print(head(Data.part))
    print("Length of moving average data:")
    print(length(Data.part[,2]))
    # set max x axis for plot by finding the last marker position
    max_x_lim = max(Data.part[,2])
    # save locations of NAs in 3rd column of data frame ("Fst.Fct")
    where.NA = which(is.na(Data.part[,3]))
    #
    # subset het data to only include that chromosome and save as separate data frame
    Het1.part = subset(het1, subset=chromosome %in% ijk)
    Het2.part = subset(het2, subset=chromosome %in% ijk)
    print(head(Het1.part))
    print("Length of moving average data, het 1:")
    print(length(Het1.part[,2]))
    print(head(Het2.part))
    print("Length of moving average data, het 2:")
    print(length(Het2.part[,2]))
    # save locations of NAs in 3rd column of data frame ("Het")
    where.NA1 = which(is.na(Het1.part[,3]))
    where.NA2 = which(is.na(Het2.part[,3]))
    #
    # if you want to export, create file name
    if( export == TRUE ){
      plotname = paste(name, ijk, sep="_")
      png(paste(plotname, "png", sep="."), width=960, height=600)
    }
    #
    # create plot layout
    layout(matrix(c(1,2), nrow = 2, byrow = TRUE), heights=c(3,1.75))
    ## PLOT FST FIRST ##
    # use min/max_lim above to define the y axis
    par(mar=c(1,4,4,2) + 0.1)
    plot(Data.part[,2], Data.part[,3],type="l",ylim=c(0,1.0), xlim=c(0,max_x_lim), lwd=3, col="red", ylab='',xlab='',xaxt = 'n', cex.lab=1, cex.axis=1)
    # add title, lines to plot
    lines(Data.part[,2], Data.part[,5], col=rgb(0.193,0.205,0.205,0.25))
    lines(Data.part[,2], Data.part[,6], col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)!=0) polygon(c(Data.part[,2][-where.NA], rev(Data.part[,2][-where.NA])), c(Data.part[,5][-where.NA], rev(Data.part[,6][-where.NA])), col=rgb(0.193,0.205,0.205,0.25))
    if(length(where.NA)==0) polygon(c(Data.part[,2], rev(Data.part[,2])), c(Data.part[,5], rev(Data.part[,6])), col=rgb(0.193,0.205,0.205,0.25))    
    abline(h=0)
    #mtext(paste(ijk),outer=TRUE,line=-2,cex=2,at=0.515)
    #mtext("Map Position (bp)",outer=TRUE,side=1,line=-2,cex=2,at=0.525)
    #mtext(expression(paste("F"[ST])),outer=TRUE,side=2,line=-2,cex=2,at=0.525)
    ## PLOT HET SECOND ##
    par(mar=c(5,4,0,2) + 0.1)
    plot(Het1.part[,2], Het1.part[,3],type="l",ylim=c(0,0.5), xlim=c(0,max_x_lim), lwd=3, col="green", ylab='',xlab='',cex.lab=1, cex.axis=1, yaxs="i")
    ## plot second data set. use min/max_lim above to define the x axis
    par(new=TRUE)
    plot(Het2.part[,2], Het2.part[,3],type="l",ylim=c(0,0.5), xlim=c(0,max_x_lim), lwd=3, col="blue", ylab='',xlab='',xaxt = "n", yaxt = "n", yaxs = "i")
    # add title, lines to plot
    abline(h=0)
    mtext(paste(ijk),outer=TRUE,line=-2,cex=2,at=0.515)
    mtext("Map Position (bp)",outer=TRUE,side=1,line=-2,cex=2,at=0.525)
    mtext(expression(paste("H"[observed])),outer=TRUE,side=2,line=-2,cex=1.5,at=0.225)
    mtext(expression(paste("F"[st])),outer=TRUE,side=2,line=-2,cex=2,at=0.65)
    legend("topright", legend=c(legend.text.het1, legend.text.het2), fill = c("green", "blue"), bty = "n", cex = 1)
    # write out the plot 
    if( export == TRUE ){
      dev.off()
    }
  }
}
  
