### PREPARE for ADDITIONAL FILTERING SCRIPTS: Remove the heading from the catalog.snps.tsv file, and add a "_" between Cat and ID in the batch.haplotypes.tsv file ###

## Arguments: (1) catalog.snps file (absolute path preferred), (2) name for new catalog.snps file (absolute path preferred), (3) haplotypes.tsv file (absolute path preferred), (4) name for new haplotypes.tsv file (absolute path preferred). Be sure to UNZIP the appropriate files first. 

import sys
snps = open(sys.argv[1], "r")
new_snps = open(sys.argv[2], "w")
haps = open(sys.argv[3], "r")
new_haps = open(sys.argv[4], "w")

#open the snps file, read all lines into new file except the first line. 
new_list = snps.readlines()[1:]
new_str = ''.join(new_list)
new_snps.write(new_str)
snps.close()
new_snps.close()

#open the haplotypes file, add in the ""_" between Cat and ID in the first line. 
first_line = haps.readline()
new_first_line = first_line.replace("Catalog ID", "Catalog_ID")
new_haps.write(new_first_line)
rest_list = haps.readlines()[1:]
rest_str = ''.join(rest_list)
new_haps.write(rest_str)
haps.close()
new_haps.close()
