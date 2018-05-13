test_df_out <- read_delim("D:/Pacific cod/DataAnalysis/PCod-Compare-repo/analyses/environmental/test_df_out.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

colnames(test_df_out) <- c("lat", "temp_125", "temp_125.5", "depth")

ggplot(test_df_out, aes(x=lat,y=temp_125)) +
  geom_point(aes(col=depth), size = 3)
