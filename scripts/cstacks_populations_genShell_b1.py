###### Generate Shell Script to Run pstacks --> populations ######

## MF 3/10/2017
## For US Cod Data

## Edited by MF 5/3/2017 for Korea PCod Data
## Edited 6/25/2017 for Korea PCod Data
## Edited 7/5/2017



## At Command Line: python cstacks_populations_genShell_7-6 ARG1 ARG2
##---- ARG1 = complete sample list file
##---- ARG2 = sample list for building cstacks catalog


#############################################################################


import sys

sampfilename = sys.argv[1]
catfilename = sys.argv[2]

newfile = open("cstacks_populations_8-24.sh", "w")

newfile.write("#!/bin/bash")

#cstacks
newfile.write("\n"+"#cstacks"+"\n")
catFile = open(catfilename, "r")

filestring = "cstacks -b 1 "
for line in catFile: 			#for each sample file listed in the cstacks catalog file
	sampID = line.strip()	
	if sampID.startswith("#"): 
		newstring = ""
	else: 		
		newstring = "-s ../stacks_b1_wgenome/" + sampID + " "
	filestring += newstring
catFile.close()
filestring += "-o ../stacks_b1_wgenome -g -p 6 2>> cstacks_out_b1_wgenome"
newfile.write(filestring)



newfile.write("\n\n")



##sstacks: run by line
newfile.write("\n"+"#sstacks"+"\n")
samplefile = open(sampfilename, "r")

for line in samplefile: 			#for each line in the barcode file
	linelist=line.strip().split()
	newstring = "sstacks -b 1 -c ../stacks_b1_wgenome/batch_1 -s ../stacks_b1_wgenome/" + linelist[0] + " -o ../stacks_b1_wgenome -p 6 2>> sstacks_out_b1_wgenome"	#creates a new -s entry for that sample input file
	newfile.write(newstring + "\n")		# appends new -s string to "filestring"
samplefile.close()



newfile.write("\n\n")



##populations
newfile.write("populations -b 1 -P ../stacks_b1_wgenome -M PopMap_combo.txt -t 36 -r 0.80 -p 3 -m 10 --genepop --fasta 2>> populations_out_b1_wgenome")



