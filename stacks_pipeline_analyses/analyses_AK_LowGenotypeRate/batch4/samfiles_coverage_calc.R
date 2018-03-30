library(readr)
PO010715_02_depth <- read_delim("D:/Pacific cod/DataAnalysis/PCod-Korea-repo/stacks_b7_wgenome/PO010715_02_depth.txt",
"\t", escape_double = FALSE, col_names = FALSE,
trim_ws = TRUE)

View(PO010715_02_depth)

sample_cat_loci = c()
sample_depths = c()
for (i in seq(1,3534522,1)) {
  catlocus = PO010715_02_depth[i,1]
  if ((catlocus %in% sample_cat_loci) == FALSE) {
    sample_cat_loci = c(sample_cat_loci, catlocus)
    depth = PO010715_02_depth[i,3]
    sample_depths = c(sample_depths, depth)
  }
}

output1 = paste("PO010715_02", length(sample_cat_loci), paste(sample_depths, collapse = ", "), sep = "\t")

outfile = "KOR_catalogloci_depths.txt"
cat(output1, file = outfile, append = TRUE)