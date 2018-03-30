### CONVERT filtered and transposed GENEPOP TO CSV FILE FOR MAF FILTERING: replace all tabs with commas ###
## MOSTLY For use with BASH_finalfilters.py ##

## argument 1 - absolute path to stacks files
## argument 2 - batch #

#MF 2/13/2017

###########################################

import sys
stacksDIRECT = sys.argv[1]
BATCH = sys.argv[2]

genepop = open(stacksDIRECT + "/batch_" + BATCH + ".CorrectedGenotypes_biallelic_genepop_TRANSPOSED.txt", "r") 
new_genepop = open(stacksDIRECT + "/batch_" + BATCH + ".CorrectedGenotypes_biallelic_genepop_TRANSPOSED.csv", "w") 

while True: 
	linestr = genepop.readline()
	new_linestr = linestr.replace(" ", ",")
	new_genepop.write(new_linestr)
	if not linestr: break
genepop.close()
new_genepop.close()
