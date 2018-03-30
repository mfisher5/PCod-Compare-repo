########### Filter Individuals for Missing Data ############

## MF 8/29/2017

## Written for post stacks filtering with write_random_snp populations model

####################################################################################

import argparse 

parser = argparse.ArgumentParser(description="filter individuals for missing data. assumes input file has individuals as column headers")

parser.add_argument("-f", "--input", help="genotype file in 2 x 2 matrix")
parser.add_argument("-og", "--output_good", help="output file containing 'good' loci")
parser.add_argument("-ob", "--output_bad", help="output file containing 'bad' loci")
parser.add_argument("-op", "--output_proportions", help="output file containing the missing data per individual")
parser.add_argument("-s", "--stacks_path", help="path to the directory containing your stacks files")
parser.add_argument("-p", "--percent", help="threshold to remove missing data: ___ percent or more loci missing genotypes")

args = parser.parse_args()

#open infile and temporary file which will contain the inverted matrix
infile = open(args.stacks_path + "/" + args.input, "r")
tempfile = open(args.stacks_path + "/temp_transposed_genotypes.txt", "w")

# transpose the matrix so that rows are individuals, columns are loci. 
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
print "calculating which individuals are missing too much data..."
propfile = open(args.stacks_path + "/" + args.output_proportions, "w")
propfile.write("individual\tn_missing\tp_missing\n")


tempfile = open(args.stacks_path + "/temp_transposed_genotypes.txt", "r")

good_individs = ""
good_count = 0
bad_individs = ""
bad_count = 0
linecount = 0
for line in tempfile:
	if linecount == 0:
		good_individs += line
		bad_individs += line
		linecount += 1
	else:
		individ = line.strip().split(" ")[0]
		genos = line.strip().split(" ")[1:]
		geno_list = []
		for genotype in genos:
			geno_list.append(genotype)
		total = len(geno_list)
		n_missing = len([i for i in geno_list if i == "0000"])
		p_missing = float(n_missing)/float(total)
		propfile.write(individ + "\t" + str(n_missing) + "\t" + str(p_missing) + "\n")
		if p_missing <= float(args.percent):
			good_individs += line
			good_count += 1
			linecount += 1
		else:
			bad_individs += line
			bad_count += 1
			linecount += 1
tempfile.close()
propfile.close()
print "Deleting temporary file..."
import os
os.remove(args.stacks_path + "/temp_transposed_genotypes.txt")

#write contents to files
print "writing individuals to files..."
outgood = open(args.stacks_path + "/" + args.output_good, "w")
outbad = open(args.stacks_path + "/" + args.output_bad, "w")
outgood.write(good_individs)
outbad.write(bad_individs)

outbad.close()
outgood.close()

print "done."

print "Total individuals processed: ", linecount - 1
print "Individuals retained: ", good_count
print "Individuals removed: ", bad_count
print "Missing data info per individual can be found in the 'proportions' output file."





