### Make 2 x 2 matrix out of 'populations' genepop file ###

## MF 8/25/2017
## Step 1 for Post-Stacks filtering with random SNP


#################################################################################

import argparse 

parser = argparse.ArgumentParser(description="convert populations genepop file into matrix for further filtering")

parser.add_argument("-g", "--genepop", help="genepop file from populations")
parser.add_argument("-o", "--output", help="output file")
parser.add_argument("-s", "--stacks_path", help="path to the directory containing your stacks files")

args = parser.parse_args()

#open files
infile = open(args.stacks_path + "/" + args.genepop, "r")
outfile = open(args.stacks_path + "/" + args.output, "w")

#write "sample" as the [1][1] in the output matrix
outfile.write("sample\t")

#print first line to skip over the title in the genepop
title = infile.readline()
print "working with genepop file: ", args.genepop
print title

#read in the comma-delimited list of loci, write out a tab-delimited list of loci
print "writing loci to output file..."
loci_str = infile.readline()
loci_list = loci_str.split(",")
loci_output = "\t".join(loci_list)

outfile.write(loci_output)

#write out each line with genotypes, removed the "," after the sample name
print "writing genotypes to output file..."
for line in infile:
	if not line.startswith("pop"):
		linelist = line.split(",")
		newstr = linelist[0] + linelist[1]
		outfile.write(newstr)
infile.close()
outfile.close()
print "done."





