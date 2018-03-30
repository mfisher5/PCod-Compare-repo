### This python script will create a list of loci from the `populations` output genepop file ###

# MF 8/30/2017

## Written for write_random_snp post stacks filtering

##################################################################################


import argparse

parser = argparse.ArgumentParser(description="create fasta file out of a 2x2 matrix with loci as column headers, individuals as rows")

parser.add_argument("-mat", "--matrix", help="filtered matrix file. Separation must be tab or comma")
parser.add_argument("-pos", "--loci_pos", help="loci position in matrix. [rows/columns]")
parser.add_argument("-cat", "--catalog", help="stacks batch.catalog.tags.tsv file")
parser.add_argument("-fasta", "--fasta_output", help="output file, in fasta format")
parser.add_argument("-p", "--path", help="path to stacks files, including matrix and catalog input files")

args = parser.parse_args()

import sys

#open the matrix file
infile = open(args.path + "/" + args.matrix, "r")


#create a list of loci from the matrix file (note - this must exclude the "sample" heading in the position [0][0])

print "Reading loci from file:"
print args.matrix

loci_list = []

if args.loci_pos == "columns":
	loci_snp_list = infile.readline().strip().split()[1:] # change to .split(" ") if separator is space
	for i in loci_snp_list:
		new_locus = i.split("_")[0]
		loci_list.append(new_locus)
	print "--"
	print "Discovered ", len(loci_list), " loci in matrix file."
	print "--"
elif args.loci_pos == "rows":
	infile.readline()
	for line in infile:
		loci_list.append(line.strip().split()[0].split("_")[0])
	print "--"
	print "Discovered ", len(loci_list), " loci in matrix file."
	print "--"
else:
	print "ERROR: unrecognized loci position type."
infile.close()



#open the UNZIPPED catalog file
catfile = open(args.path + "/" + args.catalog, "r")



#extract the sequences from the catalog file
print "Reading sequences from file:"
print args.catalog

fasta = ""
for line in catfile: 
	linelist = line.strip().split()
	if linelist[2] in loci_list: 
		newline = ">"+linelist[2]+"\n"+linelist[9]+"\n"
		fasta += newline
catfile.close()


#open a new file to write into 
print "Writing new fasta file..."
newfile = open(args.path + "/" + args.fasta_output, "w")
newfile.write(fasta)
newfile.close()

print "Done."
