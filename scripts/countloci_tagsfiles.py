### This script counts the number of loci in each ustacks tags file ###


import argparse 

parser = argparse.ArgumentParser(description="count number of consensus seqs in .tags files")

parser.add_argument("-s", "--samples", help="file with list of samples")
parser.add_argument("-d", "--directory", help="stacks directory with tags files; can be local path, don't include final '/'")
parser.add_argument("-o", "--output", help="output file name with local path")

args = parser.parse_args()

samplefile= open(args.samples, "r")
outfile = open(args.output, "w")
outfile.write("sample\tnum_loci\n")


for line in samplefile:
	sample = line.strip().split()[0]
	filename = args.directory + "/" + sample + ".tags.tsv"
	loci_count = 0
	tagsfile = open(filename, "r")
	print "Counting tags in sample " + sample
	for line in tagsfile:
		if "consensus" in line:
			loci_count += 1
	tagsfile.close()
	outfile.write("\n" + sample + "\t" + str(loci_count))

samplefile.close()

outfile.close()


