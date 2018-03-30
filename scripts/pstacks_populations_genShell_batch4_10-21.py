###### Generate Shell Script to Run pstacks --> populations ######

## MF 3/10/2017
## For US Cod Data

## Edited by MF 5/3/2017 for Korea PCod Data
## Edited 6/25/2017 for Korea PCod Data
## Edited 10/9/2017 for combined data - KOR and AK. DOES INCLUDE CSTACKS



## At Command Line: python cstacks_populations_genShell_m3_10-9.py ARG1 ARG2 ARG3
##---- ARG1 = samples for cstacks file (txt)
##---- ARG2 = Korea sample list file (population map)
##---- ARG3 = Alaska sample list file (population map)


#############################################################################


import sys



newfile = open("cstacks_populations_batch4_10-21.sh", "w")
newfile.write("#!/bin/bash\n")



#cstacks
newfile.write("\n"+"#cstacks"+"\n")
catFile = open(sys.argv[1], "r")

filestring = "cstacks -b 4 "
for line in catFile: 			#for each sample file listed in the cstacks catalog file
	sampID = line.strip()	
	if sampID.startswith("#"): 
		newstring = ""
	else: 		
		newstring = "-s ../stacks_b4_wgenome/" + sampID + " "
	filestring += newstring
catFile.close()
filestring += "-o ../stacks_b4_wgenome -g -p 6 2>> ../stacks_b4_wgenome/cstacks_out_b4_wgenome"
newfile.write(filestring)



##sstacks: run by line
newfile.write("\n"+"#sstacks"+"\n")
korsamples = open(sys.argv[2], "r")

for line in korsamples: 			#for each line in the barcode file
	linelist=line.strip().split()
	newstring = "sstacks -b 4 -c ../stacks_b4_wgenome/batch_4 -s ../stacks_b4_wgenome/" + linelist[0] + " -o ../stacks_b4_wgenome -p 6 2>> ../stacks_b4_wgenome/sstacks_out_b4_wgenome"	#creates a new -s entry for that sample input file
	newfile.write(newstring + "\n")		# appends new -s string to "filestring"
korsamples.close()

aksamples = open(sys.argv[3], "r")

for line in aksamples: 			#for each line in the barcode file
	linelist=line.strip().split()
	newstring = "sstacks -b 4 -c ../stacks_b4_wgenome/batch_4 -s ../stacks_b4_wgenome/" + linelist[0] + " -o ../stacks_b4_wgenome -p 6 2>> ../stacks_b4_wgenome/sstacks_out_b4_wgenome"	#creates a new -s entry for that sample input file
	newfile.write(newstring + "\n")		# appends new -s string to "filestring"
aksamples.close()


newfile.write("\n\n")



##populations
newfile.write("populations -b 4 -P ../stacks_b4_wgenome -M PopMap_combo_b4.txt -t 36 -r 0.5 -p 3 -m 10 --genepop --fasta 2>> ../stacks_b4_wgenome/populations_out_b4_wgenome")



