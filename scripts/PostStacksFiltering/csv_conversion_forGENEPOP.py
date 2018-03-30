#### This script converts the FINAL filtered .csv document (corrected genotypes, MAF filtered, low coverage loci filtered, low coverage individuals filtered) into a genepop file format ####

## MF 8/28/2017

## Based on "genepop conversion for R" script MF 12/13/2016
## Sample names should be column headers, Loci should be rows


#################################################################################

import argparse 

parser = argparse.ArgumentParser(description="convert matrix file (samples as column headers) to a genepop file.")

parser.add_argument("-f", "--input", help="2 x 2 matrix file with samples as column headers, one row per locus.")
parser.add_argument("-o", "--output", help="output file in genepop format")
parser.add_argument("-m", "--popmap", help="population map file")
parser.add_argument("-pm", "--path_popmap", help="path to population map file")
parser.add_argument("-ps", "--path_stacks", help="path to stacks files, including input file")
parser.add_argument("-t", "--title", help="title for your genepop file")
parser.add_argument("-split", "--split_by", help="split contents of each line by... ['tab' / 'space']")


args = parser.parse_args()

################## This first part of the script uses the popmap file to generate the column indices needed for the rest of the script. #####################
## Generate list of populations, and number of samples per population
import sys
popmap = open(args.path_popmap + "/" + args.popmap, "r")

#--- use of collections module keeps the Dictionary keys in the order in which they are entered - so that the order matches `populations` output files. THIS IS IMPORTANT
import collections
PopDict = collections.OrderedDict()
PopList = []
for line in popmap: 
	linelist = line.strip().split()
	newpop = linelist[1]
	if newpop not in PopDict: 
		PopDict[newpop] = 1
		PopList.append(newpop)
	elif newpop in PopDict: 
		count = PopDict[newpop]
		count += 1
		PopDict[newpop] = count

#--- print message
length = str(len(PopList))
print "You have " + length + " populations."
print "These are your populations, with the number of samples in each:" 
print PopDict



#--- create column counts
column = 1
column_indices = {}
for pop in PopList: 
	n_samples = PopDict[pop]
	new_column = column + n_samples
	column_indices[pop] = [column, new_column]
	column = new_column




print ("creating script for part 2...")


################## This second part of the script writes a second script, which will be called with 'subprocess,' to write the genepop file ####################

part2 = open("csv_conversion_forGENEPOP_p2.py", "w")

part2.write("### This is part two of the script that will convert a 2x2 matrix file into a genepop file. It is automatically generated with part 1 #####\n\n\n")
part2.write("infile = open('" + args.path_stacks + "/" + args.input + "', 'r')\n")
part2.write("genepop = open('" + args.path_stacks +  "/" + args.output + "', 'w')\n")

# create a title for the genepop file
part2.write("genepop.write('" + args.title + r"\r\n')" + "\n") #An "r" before the quotes allows you to write "\n" without creating a new line
 

part2.write("print 'transposing genotypes matrix...'\n")
# transpose the matrix file so that the loci are along the top row, and the individual names are in the first column
part2.write("data_matrix = []\n")
if args.split_by == "space":
	part2.write("for line in infile:\n\ttmp_line = ''\n\ttmp_line += line\n\tdata_matrix.append(tmp_line.split(' '))\n")
elif args.split_by == "tab":
	part2.write("for line in infile:\n\ttmp_line = ''\n\ttmp_line += line\n\tdata_matrix.append(tmp_line.split" + r"('\t'))" + "\n")
elif args.split_by == "comma":
	part2.write("for line in infile:\n\ttmp_line = ''\n\ttmp_line += line\n\tdata_matrix.append(tmp_line.split" + r"(','))" + "\n")
else:
	print "ERROR: incorrect split argument"
part2.write("infile.close()\n\n")
part2.write("transposed = zip(*data_matrix)\n\n")


part2.write("print 'writing loci into genepop file...'\n")
#create loci list and write it to the genepop file.
part2.write("locilist = transposed[0]\n")
part2.write("LociIndex = range(0, len(locilist))\n")
part2.write("for i in LociIndex:\n\tif transposed[0][i] != 'sample':\n\t\tgenepop.write(transposed[0][i] + " + r"'\r\n')" + "\n\n")


# generate a string that splits by columns for each population and write to script
n_pops = len(PopList)
pop_count = 1
for pop in PopList: 
	if pop_count < n_pops:
		cols = column_indices[pop]
		newstr = pop + " = transposed[" + str(cols[0]) + ":" + str(cols[1]) + "]"
		part2.write(newstr + "\n")
		pop_count += 1
	elif pop_count == n_pops:
		cols = column_indices[pop]
		last_column = cols[1] - 1
		newstr = pop + " = transposed[" + str(cols[0]) + ":" + str(last_column) + "]"
		part2.write(newstr + "\n")
		part2.write("last_line = list(transposed[" + str(last_column) + "])\n")
		part2.write("seq = range(0, len(last_line))\n")
		part2.write("for i in seq:\n\tlast_line[i] = last_line[i].strip(" + r"'\r\n')" + "\n")
part2.write("\n\n")

part2.write("print 'writing genotypes into genepop file by population...'\n")
# write the genotypes to a genepop file by population # 

for pop in PopList:
	part2.write("genepop.write('Pop' + " + r"'\r\n')" + "\n")
	part2.write("for line in " + pop + ":\n\t")
	part2.write("linestr = " + r"'\t'" + ".join(line[1:])\n\t")
	part2.write("newline = line[0] + " + r"',\t' + linestr + '\r\n'" + "\n")
	part2.write("\tgenepop.write(newline)\n\n")

part2.write("linestr = '\t'.join(last_line[1:])\n")
part2.write("newline = last_line[0] + " + r"',\t'" + " + linestr + " +  r"'\r\n'" + "\n")
part2.write("genepop.write(newline)\n\n")


part2.write("genepop.close()\n\n")
part2.write("print 'done.'")

part2.close()


print "calling script for part 2..."

import subprocess

subprocess.call(['python', 'csv_conversion_forGENEPOP_p2.py'])


