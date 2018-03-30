############# This script parses genepop output from the .INF file (Option 5 >> Suboption 1) ##################

## MF 8/29/2017
## Written for MAF filtering 

###################################################################################

import argparse

parser = argparse.ArgumentParser(description="parse .INF file from genepop to get allele frequences by population by locus")

parser.add_argument("-f", "--input", help=".INF file")
parser.add_argument("-p", "--path_stacks", help="path to stacks files, including .INF file")
parser.add_argument("-pop", "--pops", help="text file with population names, one on each line. make sure this is in the same order as the populations listed in the .INF output file! (aka the order of your populations in your pop map file)")
parser.add_argument("-b", "--batch", help="stacks batch number")

args = parser.parse_args()


# open files 
infile = open(args.path_stacks + "/" + args.input, "r")
outfile = open(args.path_stacks + "/batch_" + args.batch + "_parseINF.txt", "w")

#### preparing for the parsed file ####


# skip down to where the tables of allele frequencies start
for line in infile:
	if line.startswith("Tables of allelic frequencies for each locus"):
		break

# for each locus, save the number of alleles IF that is the greatest number of alleles at any locus 
# this will help generate the header of the file. 
max_n_alleles = 0
count = 0

for line in infile:
	if line.startswith(" Locus") and count == 0:
		count = 1
	elif count == 1 or count == 2 or count == 3:
		count += 1
	elif count == 4: 
		linelist = line.strip().split(" ")
		n_alleles = len([i for i in linelist if i != ""])
		if n_alleles > max_n_alleles:
			max_n_alleles = n_alleles
		count = 0
infile.close()
print "You have a maximum of ", max_n_alleles, " alleles at one or more of your loci."

#create the header for the frequency output file
header = "locus\t"

pops_file = open(args.path_stacks + "/" + args.pops, "r")
pops_list = []
for line in pops_file:
	pops_list.append(line.strip().split()[0])
print len(pops_list), " populations detected."
pops_file.close()

for pop in pops_list:
	for i in range(1, max_n_alleles+1):
		header += pop + "_Allele" + str(i) + "\t"

outfile.write(header + "\n")

print "Heading written to output file. Now parsing allele frequencies..."

#### parsing the allele frequencies ####


# re-open infile
infile = open(args.path_stacks + "/" + args.input, "r")

# parse and rewrite allele frequencies

for line in infile:
	if line.startswith("Tables of allelic frequencies for each locus"):
		break

n_pops = len(pops_list)
last_pop = n_pops + 4
loci_count = 0

for line in infile:
	if line.startswith(" Locus") and count == 0:
		linelist = line.strip().split(" ")
		new_locus = linelist[1]
		newstr = new_locus
		count = 1
		loci_count += 1
	elif count == 1 or count == 2 or count == 3:
		count += 1
	elif count == 4: 
		linelist = line.strip().split(" ")
		n_alleles = len([i for i in linelist if i != ""])
		count += 1
	elif count >= 5 and count < last_pop:
		linelist = line.strip().split(" ")
		new_pop = linelist[0]
		allele_info = [i for i in linelist[1:] if i != '']
		for i in range(0,n_alleles):
    			newstr += "\t" + allele_info[i]
		if n_alleles < max_n_alleles:
    			alleles_to_add = max_n_alleles - n_alleles
    			for x in range(0,alleles_to_add):
        			newstr += "\t"
		count += 1
	elif count == last_pop:
		linelist = line.strip().split(" ")
		new_pop = linelist[0]
		allele_info = [i for i in linelist[1:] if i != ""]
		for i in range(0,n_alleles):
    			newstr += "\t" + allele_info[i]
		if n_alleles < max_n_alleles:
    			alleles_to_add = max_n_alleles - n_alleles
    			for x in range(0,alleles_to_add):
        			newstr += "\t"
		outfile.write(newstr + "\n")
		count = 0
infile.close()
outfile.close()

print "Done."
print "Parsed allele frequencies in ", loci_count, " loci."
	

