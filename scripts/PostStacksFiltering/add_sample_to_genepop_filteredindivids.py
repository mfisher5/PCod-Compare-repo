#### This script adds in the word "sample" to your column headers for the genepop file, filtered by Marine's scripts###

### Mostly for use with BASH_finalfilters. ASSUMES THIS FILE NAMING SCHEME ###

## argument 1 - Absolute path to stacks files (or wherever the filtered genepop file is), without final "/"
## argument 2 - batch #
## argument 3 - filename

##########################################################################

import sys
stacksDIRECT = sys.argv[1]
BATCH = sys.argv[2]
FILENAME = sys.argv[3]

genepop = open(stacksDIRECT + "/" + sys.argv[3], "r") 
new_genepop = open(stacksDIRECT + "/batch_" + BATCH + ".CorrectedGenotypes_biallelic_edit.genepop", "w") 

first_line = genepop.readline()
if first_line.startswith("sample"): 
	print "original genepop file ok"
	temp_list = genepop.readlines()
	temp_str = ''.join(temp_list)
	new_genepop.write(temp_str)
else: 
	print "adding 'sample' as the first column header"
	new_first_line = "sample" + first_line
	new_genepop.write(new_first_line)
	rest_list = genepop.readlines()
	rest_str = ''.join(rest_list)
	new_genepop.write(rest_str)

genepop.close()
new_genepop.close()
