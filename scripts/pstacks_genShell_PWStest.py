###### Generate Shell Script to Run pstacks --> populations ######

## MF 3/10/2017 Edited 5/25/2017
## For US Cod Data



## At Command Line: python cstacks_populations_genShell_3-8 ARG1 ARG2
##---- ARG1 = complete sample list file (barcode file)


#############################################################################


import sys

sampfilename = sys.argv[1]

newfile = open("pstacks_PWStest.sh", "w")
newfile.write("#!/bin/bash")

#pstacks

newfile.write("\n"+"#pstacks"+"\n")
samplefile = open(sampfilename, "r")
ID_int = 707
for line in samplefile: 			#for each line in the barcode file	
	sampID = line.strip().split()[1]	
	ustacks_code = "pstacks -t sam -f ../stacks_b1_wgenome/" + sampID + ".sam -o ../stacks_b1_wgenome -i " + str(ID_int) + " -m 10 -p 6 --model_type bounded 2>> pstacks_out_b1_wgenome_PWStest" + "\n"
	newfile.write(ustacks_code)	#append this new line of code to the output file
	ID_int += 1
samplefile.close()

newfile.close()


