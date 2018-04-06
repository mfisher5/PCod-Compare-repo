################ Filter Sliding Window Analysis by Marker Density #################
#
# Use this script to remove data from sliding window analysis where moving average was...
#  calculated with fewer than `x` markers. 
#
# MF 3/24/2018
#
########################################################################################


# Function --------------------------------------------------------
filter_marker_density = function(swa_output = swa_output, swa_input = swa_input, window_size = 250000, divisions = 150, cutoff = 2, outfile = "output_kernel_smoothing_1e+5_bootstraps_sigma_250000_div15_FILTERED.txt"){
  print("Warning: this function assumes that all positions are unique")
  window_size=window_size
  divisions=divisions
  Num_marker_data = nrow(swa_input)
  Num_chromosome = unique(swa_input$chromosome)
  chromosomes_all= rep(Num_chromosome,each=divisions) # `each` argument should == divisions
  Positions_all=c()
  markers_per_window=c()
  ## CHARLIE'S CODE: for loop to calculate markers per window on each chromosome
  for (ijx in Num_chromosome){
    Data.chrom = subset(swa_input, swa_input$chromosome==ijx)
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
    }
  }
  ## save output as a matrix
  output=c()   ##be sure to clear output from previous runs
  output=data.frame(chromosomes_all,Positions_all,markers_per_window)
  colnames(output)=c("chrom","position","num_markers")
  
  marker_density = output
  cutoff = cutoff
  
  ## remove data when num_markers is less than cutoff
  marker_density$num_markers[marker_density$num_markers < cutoff] <- NA
  marker_density$position <- as.character(marker_density$position)
  pos_to_remove <- marker_density[is.na(marker_density$num_markers),2]
  print("Number of data points to remove based on marker density:")
  print(length(pos_to_remove))
  swa_output$position <- as.character(swa_output$position)
  count = 0
  swa_output_filtered <- swa_output
  for(i in swa_output$position){
    newpos <- i
    if(newpos %in% pos_to_remove){
      count = count + 1
      row_index <- which(swa_output$position == i)
      swa_output_filtered[row_index,3] <- NA
    }
  }
  print("Number of positions actually removed from analysis:")
  print(count)
  write.table(swa_output_filtered, outfile, sep ="\t", row.names = FALSE, quote=FALSE)
}




# Read in Data ------------------------------------------------------------
swa_output = read.table("EastvWest/batch_8_SWA_eastwest_globalFst_filtered_output_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150.txt", header = TRUE, sep = "\t")
swa_input = read.table("EastvWest/batch_8_SWA_input_eastwest_globalFst_filtered_sorted.txt", header = TRUE, sep = "\t")



swa_output_west = read.table("West/batch_8_final_filtered_west_2reg_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150.txt", header = TRUE, sep = "\t")
swa_input_west = read.table("batch_8_SWA_input_west_sorted.txt", header = TRUE, sep = "\t")

swa_output_east = read.table("East/batch_8_final_filtered_east_globalFST_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150.txt", header = TRUE, sep = "\t")
swa_input_east = read.table("East/batch_8_SWA_input_east_sorted.txt", header = TRUE, sep = "\t")




# Run Function ------------------------------------------------------------
filter_marker_density(swa_output = swa_output, swa_input = swa_input, window_size = 250000, divisions = 150, cutoff = 2, outfile = "EastvWest/batch_8_SWA_eastwest_globalFst_filtered_output_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150_FILTERED.txt")


filter_marker_density(swa_output = swa_output_west, swa_input = swa_input_west, window_size = 250000, divisions = 150, cutoff = 2, outfile = "West/batch_8_final_filtered_west_2reg_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150_FILTERED.txt")



filter_marker_density(swa_output = swa_output_east, swa_input = swa_input_east, window_size = 250000, divisions = 150, cutoff = 2, outfile = "batch_8_final_filtered_east_globalFST_kernel_smoothing_1e+05_bootstraps_sigma_250000_div150_FILTERED.txt")


