#### This script will take a text file containing the R output from BAYESCAN's function plot_bayescan() and find the correct outlier loci's IDs ####

## MF 10/6/2017
## Note that input file should be a text file of R output in the console copied directly into a text file 

#######################################################################################################################
import argparse 

parser = argparse.ArgumentParser(description="Match bayescan outlier loci IDs to the actual stacks IDs (if PGD spider was used for file conversion).")

parser.add_argument("-i", "--input", help="fst text file output from bayescan")
parser.add_argument("-gen", "--genepop", help="the genepop file used in PGD spyder to create BAYESCAN input file")
parser.add_argument("-o", "--output", help="output text file")
parser.add_argument("-s", "--separator", help="separator used in genepop file [comma/newline]")


args = parser.parse_args()

infile = open(args.input, "r")
genepop = open(args.genepop, "r")
outfile = open(args.output, "w")


## put stacks loci in a list so indices will correspond to bayescan locus IDs 
print "indexing stacks loci..."
stacks_loci = []

genepop.readline() #header of genepop

line = genepop.readline()

if args.separator == "comma":
	stacks_loci = line.strip().split(",")
elif args.separator == "newline":
	while not line.startswith("Pop"):
		stacks_loci.append(line.strip())
		line = genepop.readline()
else:
	print "ERROR: unacceptable separator argument."
genepop.close()

print "You have ", len(stacks_loci), " loci."


## create header of output file with bayescan header 
### will get rid of the weird spacing in bayescan and have first column listed as "locus"
headerlist = infile.readline().strip().split()
header = "locus " + " ".join(headerlist) + "\n"
outfile.write(header)


## read over lines from bayescan with locus names from stacks
print "copying over BAYESCAN output.."

line = infile.readline()
count = 0
while line:
	linelist = line.strip().split()[1:]
	locus = stacks_loci[count]
	newline = locus + " " + " ".join(linelist) + "\n"
	outfile.write(newline)
	line = infile.readline()
	count += 1

print "Copied over ", count, " loci."
 
infile.close()
outfile.close()