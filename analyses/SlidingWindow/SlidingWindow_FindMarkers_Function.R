############# Find Marker Density & Loci Names per Window, Sliding Window Analysis #########
#
# Adopted from Charlie Waters' script "just_plot_updated_changed_plots"
# To be run AFTER SlidingWindow_MF.R script
# Will not produce plots, only text file
#
# MF edited for Pacific cod 2/26/2018
#
#########################################################################


find_markers_in_window = function(marker_data = marker_data, window_size= 250000, divisions = 150, output = "Num_Name_Loci_Per_Window.txt"){
  window_size=window_size
  divisions=divisions
  Num_marker_data = nrow(marker_data)
  Num_chromosome = unique(marker_data$chromosome)
  chromosomes_all= rep(Num_chromosome,each=divisions) # `each` argument should == divisions
  Positions_all=c()
  markers_per_window=c()
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
      markers_per_window = c(markers_per_window,num_markers)
      window_subset <- Data.chrom[MA_window,]
      name_markers <- c(window_subset$Locus)
      name_markers_line <- paste(name_markers, collapse=",")
      newline <- paste(ijx, Positions_LG[x], num_markers, name_markers_line, sep="\t")
      write(newline, file=output, append=TRUE)
    }
  }
}

