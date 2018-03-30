### Filter loci for missing data ###

## MF 8/25/2017
## Step 2 for Post-Stacks filtering with random SNP


#################################################################################

import argparse 

parser = argparse.ArgumentParser(description="filter loci for missing data")

parser.add_argument("-f", "--input", help="genotype file in 2 x 2 matrix")
parser.add_argument("-og", "--output_good", help="output file containing 'good' loci")
parser.add_argument("-ob", "--output_bad", help="output file containing 'bad' loci")
parser.add_argument("-op", "--output_proportions", help="output file containing the missing data per locus")
parser.add_argument("-s", "--stacks_path", help="path to the directory containing your stacks files")
parser.add_argument("-p", "--percent", help="threshold to remove missing data: ___ percent or more missing genotypes")

args = parser.parse_args()

#open infile and temporary file which will contain the inverted matrix
infile = open(args.stacks_path + "/" + args.input, "r")
tempfile = open(args.stacks_path + "/temp_transposed_genotypes.txt", "w")

# Transpose the input matrix and write out to a temporary file. ## Tranpose code taken from Dan Drinan's `transpose.py` script. ##
print "Transposing matrix..."
header = True
matrix_of_data = []

for line in infile:
    tmp_line = ''

    if header:
        tmp_line = ''
        header = False

    tmp_line += line

    matrix_of_data.append(tmp_line.split())

infile.close()
        
transposed = zip(*matrix_of_data)

for line in transposed:
	tmp_output = str(line).replace('[', '').replace('(', '').replace("'", '').replace(')', '').replace(',', '').replace(']', '') + '\n'
	tempfile.write(tmp_output)

tempfile.close()


#read in the contents of the temporary file. calculate missing data for each locus, and save line to a 'good' or 'bad' string
print "calculating 'good' and 'bad' loci..."
propfile = open(args.stacks_path + "/" + args.output_proportions, "w")
propfile.write("locus\tn_missing\tp_missing\n")

tempfile = open(args.stacks_path + "/temp_transposed_genotypes.txt", "r")

good_loci = ""
good_count = 0
bad_loci = ""
bad_count = 0
linecount = 0
for line in tempfile:
	if linecount == 0:
		good_loci += line
		bad_loci += line
		linecount += 1
	else:
		locus = line.strip().split(" ")[0]
		genos = line.strip().split(" ")[1:]
		geno_list = []
		for genotype in genos:
			geno_list.append(genotype)
		total = len(geno_list)
		n_missing = len([i for i in geno_list if i == "0000"])
		p_missing = float(n_missing)/float(total)
		propfile.write(locus + "\t" + str(n_missing) + "\t" + str(p_missing) + "\n")
		if p_missing <= float(args.percent):
			good_loci += line
			good_count += 1
			linecount += 1
		else:
			bad_loci += line
			bad_count += 1
			linecount += 1
tempfile.close()
propfile.close()

#write contents to files
print "writing 'good' and 'bad' loci to files..."
outgood = open(args.stacks_path + "/" + args.output_good, "w")
outbad = open(args.stacks_path + "/" + args.output_bad, "w")
outgood.write(good_loci)
outbad.write(bad_loci)

outbad.close()
outgood.close()

print "done."

print "Total loci: ", linecount - 1
print "Loci retained: ", good_count
print "Loci removed: ", bad_count
print "Missing data info per locus can be found in the 'proportions' output file."

