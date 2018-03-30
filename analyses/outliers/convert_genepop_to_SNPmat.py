############ PRODUCES THE INPUT FILE / LISTS REQUIRED FOR OUTFLANK ##########

## This script will produce the required input for OutFlank's function MakeDiploidFSTMat(), which creates the required input file format for OutFLANK
## Input should be a genepop file with a header line, and a population map file (can include more individuals than are in the genepop file)
## Use this with the R OutFLANK script for KorPCod here: 
## MF 10/5/2017

#################################################################################


import argparse 

parser = argparse.ArgumentParser(description="produce SNPmat file, and files containing loci / population lists for OutFLANK outlier analysis.")

parser.add_argument("-i", "--input", help="genepop file that you want to run through OutFLANK")
parser.add_argument("-p", "--popmap", help="population map from stacks (each line has sample - tab - population")
parser.add_argument("-o", "--output", help="bash shell script file name. must have file extension .sh")
parser.add_argument("-ol", "--outLocusNames", help="text file with the name of each locus on each line, to be read into R")
parser.add_argument("-op", "--outPopNames", help="text file with the name of each sample's population on each line, to be read into R")


args = parser.parse_args()


import sys
infile = open(args.input, "r")
popmap = open(args.popmap, "r")
outfile = open(args.output, "w")
locusfile = open(args.outLocusNames, "w")
popfile = open(args.outPopNames, "w")


########### Make SNPmat file and get locusfile output #########
header = infile.readline()
print header


loci_str = ""
line = infile.readline()
while not line.startswith("Pop") and not line.startswith("pop"):
	loci_str += '"' + line.strip() + '"' + "\n"
	line = infile.readline()

locusfile.write(loci_str.strip())

genotypes = ""
for line in infile:
	if not line.startswith("Pop") and not line.startswith("pop"):
		tmp_genotypes_list = line.strip().split()[1:]
		tmp_genotypes = "\t".join(tmp_genotypes_list)
		newline = tmp_genotypes.replace('0000', '9').replace('0101', '0').replace('0202', '2').replace('0102', '1') + "\n"
		genotypes += newline
infile.close()

outfile.write(genotypes)
outfile.close()

print "Done creating SNPmat file."



####### Make Popfile output #########

pop_dict = {}
for line in popmap:
	pop_dict[line.strip().split()[0]] = line.strip().split()[1]
popmap.close()
	
infile = open(args.input, "r")
infile.readline() # header


pops_str = ""
for line in infile:
	linelist = line.strip().split()
	if len(linelist) > 1:
		sample = linelist[0].strip(",")
		pops_str += '"' + pop_dict[sample] + '"\n'
infile.close()

popfile.write(pops_str.strip())
popfile.close()



