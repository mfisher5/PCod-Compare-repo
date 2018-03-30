### This python script will create a list of loci from the `populations` output genepop file ###

# MF 8/30/2017

## Written for write_random_snp post stacks filtering

##################################################################################


import argparse

parser = argparse.ArgumentParser(description="create fasta file out of a 2x2 matrix with loci as column headers, individuals as rows")

parser.add_argument("-i", "--input", help="input, filtered genepop file.")
parser.add_argument("-o", "--output", help="output file, in fasta format")
parser.add_argument("-c", "--catalog", help="stacks batch.catalog.tags.tsv file")
parser.add_argument("-s", "--gen_sep", help="how are loci in the genepop file separated? [comma/newline]")
parser.add_argument("-p", "--path", help="path to stacks files, including matrix and catalog input files")

args = parser.parse_args()

import sys

#open the matrix file
infile = open(args.path + "/" + args.input, "r")
header = infile.readline()


#create a list of loci from the matrix file (note - this must exclude the "sample" heading in the position [0][0])

print "Reading loci from genepop:"
print args.input

loci_list = []

if args.gen_sep == "comma":
	loci_snp_list = infile.readline().strip().split()
	for i in loci_snp_list:
		new_locus = i.split("_")[0]
		loci_list.append(new_locus)
	print "--"
	print "Discovered ", len(loci_list), " loci in matrix file."
	print "--"
elif args.gen_sep == "newline":
	line = infile.readline()
	while not line.startswith("Pop"):
		loci_list.append(line.strip().split()[0].split("_")[0])
		line = infile.readline()
	print "--"
	print "Discovered ", len(loci_list), " loci in matrix file."
	print "--"
else:
	print "ERROR: unrecognized locus separation in genepop file."
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
newfile = open(args.path + "/" + args.output, "w")
newfile.write(fasta)
newfile.close()

print "Done."
