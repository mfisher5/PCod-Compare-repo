###### Generate Shell Script to Run pstacks --> populations ######

## MF 3/10/2017
## For US Cod Data

## Edited by MF 5/3/2017 for Korea PCod Data
## Edited 6/25/2017 for Korea PCod Data
## Edited 10/9/2017 for combined data 



## At Command Line: python cstacks_populations_genShell_7-6 ARG1
##---- ARG1 = Alaska sample list file (population map)


#############################################################################


import sys

newfile = open("pstacks_batch4_10-20_kor.sh", "w")
newfile.write("#!/bin/bash\n")

#pstacks WITH m OF 3

newfile.write("\n"+"#pstacks"+"\n")


aksamples = open(sys.argv[1], "r")
ID_int = 1
for line in aksamples: 			#for each line in the barcode file	
	sampID = line.strip().split()[0]
	if ID_int < 10:
		pstacks_code = "pstacks -t sam -f ../stacks_b4_wgenome/" + sampID + "_92bp.sam -o ../stacks_b4_wgenome -i 00" + str(ID_int) + " -m 3 -p 6 --model_type bounded --bound_high 0.05 2>> ../stacks_b4_wgenome/pstacks_out_b4_wgenome" + "\n"
	elif ID_int >= 10 and ID_int < 100:
		pstacks_code = "pstacks -t sam -f ../stacks_b4_wgenome/" + sampID + "_92bp.sam -o ../stacks_b4_wgenome -i 0" + str(ID_int) + " -m 3 -p 6 --model_type bounded --bound_high 0.05 2>> ../stacks_b4_wgenome/pstacks_out_b4_wgenome" + "\n"
	elif ID_int >= 100:
		pstacks_code = "pstacks -t sam -f ../stacks_b4_wgenome/" + sampID + "_92bp.sam -o ../stacks_b4_wgenome -i " + str(ID_int) + " -m 3 -p 6 --model_type bounded --bound_high 0.05 2>> ../stacks_b4_wgenome/pstacks_out_b4_wgenome" + "\n"
	newfile.write(pstacks_code)	#append this new line of code to the output file
	ID_int += 1
aksamples.close()


newfile.close()

