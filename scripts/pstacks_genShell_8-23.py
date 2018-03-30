###### Generate Shell Script to Run pstacks --> populations ######

## MF 3/10/2017
## For US Cod Data

## Edited by MF 5/3/2017 for Korea PCod Data
## Edited 6/25/2017 for Korea PCod Data, 8/21/2017 to include batch number as an argument

## At Command Line: python cstacks_populations_genShell_7-6 ARG1 ARG2
##---- ARG1 = complete sample list file, can be population map
##---- ARG2 = batch number


#############################################################################


import sys

sampfilename = sys.argv[1]
batch = sys.argv[2]

newfile = open("pstacks_b" + batch + ".sh", "w")

#pstacks

newfile.write("\n"+"#pstacks"+"\n")
samplefile = open(sampfilename, "r")
ID_int = 330
for line in samplefile: 			#for each line in the barcode file	
	sampID = line.strip().split()[0]	
	ustacks_code = "pstacks -t sam -f ../stacks_b" + batch + "_wgenome/" + sampID + ".sam -o ../stacks_b" + batch + "_wgenome -i " + str(ID_int) + " -m 10 -p 6 --model_type bounded 2>> ../stacks_b" + batch + "_wgenome/pstacks_out_b" + batch + "_wgenome" + "\n"
	newfile.write(ustacks_code)	#append this new line of code to the output file
	ID_int += 1
samplefile.close()

newfile.close()

