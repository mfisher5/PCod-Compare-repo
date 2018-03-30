####### This script will filter for MAF based on the parsed output file #########

## MF 8/29/2017

## Written for write_random_snp post stacks filtering

##################################################################################


import argparse

parser = argparse.ArgumentParser(description="filter for MAF using parsed .INF genepop file.")

parser.add_argument("-inf", "--INFinput", help="parsed .INF file")
parser.add_argument("-mat", "--MATRIXinput", help="half filtered matrix file that was used to make genepop file")
parser.add_argument("-og", "--output_good", help="output file, good loci")
parser.add_argument("-ob", "--output_bad", help="output file, bad loci")
parser.add_argument("-ofg", "--outputfreqs_good", help="output file, frequencies of good loci")
parser.add_argument("-ofb", "--outputfreqs_bad", help="output file, frequencies of bad loci")
parser.add_argument("-p", "--path_stacks", help="path to stacks files, including .INF file")
parser.add_argument("-a", "--alleles", help="maximum number of alleles; this should have been printed out during the parsing script.")
parser.add_argument("-b", "--batch", help="stacks batch number")

args = parser.parse_args()


# open files 
infile = open(args.path_stacks + "/" + args.INFinput, "r")
goodfreqs = open(args.path_stacks + "/" + args.outputfreqs_good, "w")
badfreqs = open(args.path_stacks + "/" + args.outputfreqs_bad, "w")

# write header lines
header = infile.readline()
goodfreqs.write(header)
badfreqs.write(header)
header_list = header.strip().split()


# identify loci with MAF < 0.05 in all populations using INF file
loci_to_keep = []
for line in infile:
	linelist = line.strip().split("\t")
	locus = linelist[0] #save locus ID
	freqs = linelist[1:] #create list of allele frequencies
	if int(args.alleles) == 2:
		allele1_freqs = freqs[0::2] #create list of allele 1 frequencies by saving every other freq in original list
		allele2_freqs = freqs[1::2] #create list of allele 2 frequencies by saving every other freq in original list, offset 1
		alleles_list = [allele1_freqs, allele2_freqs] #create list of two lists
	elif int(args.alleles0) == 3:
		allele1_freqs = freqs[0::3]
		allele2_freqs = freqs[1::3]
		allele3_freqs = freqs[2::3]
		alleles_list = [allele1_freqs, allele2_freqs, allele3_freqs]
	elif int(args.alleles) == 4:
		allele1_freqs = freqs[0::4]
		allele2_freqs = freqs[1::4]
		allele3_freqs = freqs[2::4]
		allele4_freqs = freqs[3::4]
		alleles_list = [allele1_freqs, allele2_freqs, allele3_freqs, allele4_freqs]
	count_good_alleles = 0 #reset good count to zero
	n_alleles = 0 #reset total allele count to zero
	# premise of for loop: to retain locus, must have at least one allele1 and one allele2 between 0.05 and 0.95. 
	# loop over allele 1 and allele 2 list separately
	for freq_list in alleles_list:
		# if the allele is genotyped in at least one population...
		if len([freq for freq in freq_list if freq != '-']) > 0:
			n_alleles += 1 #add to count of alleles
			genotyped_alleles = [i for i in freq_list if i != '-'] #save freqs if that allele was found in the population
			# if any allele frequencies are greater than 0.05 and less than 0.95...
			if any(float(freq) > 0.05 and float(freq) < 0.95 for freq in genotyped_alleles):
				count_good_alleles += 1 #add one to "good allele" count
	if count_good_alleles == n_alleles:
		goodfreqs.write(line)
		loci_to_keep.append(locus)
	else:
		badfreqs.write(line)

infile.close()
goodfreqs.close()
badfreqs.close()


# filter out loci based on frequencies
infile2 = open(args.path_stacks + "/" + args.MATRIXinput, "r")
goodfile = open(args.path_stacks + "/" + args.output_good, "w")
badfile = open(args.path_stacks + "/" + args.output_bad, "w")

header = infile2.readline()
goodfile.write(header)
badfile.write(header)

loci_to_good = 0
loci_to_bad = 0
for line in infile2:
	locus = line.strip().split()[0]
	if locus in loci_to_keep:
		goodfile.write(line)
		loci_to_good += 1
	else:
		badfile.write(line)
		loci_to_bad += 1

infile.close()
goodfile.close()
badfile.close()

print loci_to_good, " loci written to filtered output file."
print "Filtered out ", loci_to_bad, " loci."
	
		


