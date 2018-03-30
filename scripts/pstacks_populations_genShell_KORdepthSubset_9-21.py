###### Generate Shell Script to Run pstacks --> populations ######

## MF 3/10/2017
## For US Cod Data

## Edited by MF 5/3/2017 for Korea PCod Data
## Edited 6/25/2017 for Korea PCod Data
## Edited 7/5/2017



## At Command Line: python cstacks_populations_genShell_7-6 ARG1 ARG2
##---- ARG1 = complete sample list file


#############################################################################


import sys

sampfilename = sys.argv[1]

newfile = open("pstacks_populations_KORsubsetDepth_9-20.sh", "w")
newfile.write("#!/bin/bash\n")

#pstacks WITH m OF 5

newfile.write("\n"+"#pstacks"+"\n")
samplefile = open(sampfilename, "r")
ID_int = 1000
for line in samplefile: 			#for each line in the barcode file	
	sampID = line.strip().split()[0]	
	pstacks_code = "pstacks -t sam -f ../stacks_b2_wgenome/" + sampID + "_subset.sam -o ../stacks_b2_wgenome -i " + str(ID_int) + " -m 3 -p 6 --model_type bounded 2>> ../stacks_b2_wgenome/pstacks_out_b2_wgenome" + "\n"
	newfile.write(pstacks_code)	#append this new line of code to the output file
	ID_int += 1
samplefile.close()

newfile.write("\n\n")


##sstacks: run by line
newfile.write("\n"+"#sstacks"+"\n")
samplefile = open(sampfilename, "r")

for line in samplefile: 			#for each line in the barcode file
	linelist=line.strip().split()
	newstring = "sstacks -b 2 -c ../stacks_b1_wgenome/batch_1 -s ../stacks_b2_wgenome/" + linelist[0] + "_subset -o ../stacks_b2_wgenome -p 6 2>> sstacks_out_b2_wgenome"	#creates a new -s entry for that sample input file
	newfile.write(newstring + "\n")		# appends new -s string to "filestring"
samplefile.close()



newfile.write("\n\n")



##populations
newfile.write("populations -b 2 -P ../stacks_b2_wgenome -M PopMap_KOR_subset.txt -t 36 -r 0.8 -p 3 -m 10 --genepop --fasta 2>> populations_out_b2_wgenome")



