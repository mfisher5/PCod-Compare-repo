### This file generates Eleni's MAF filtering file ###

## MF 2/24/2017

## Arguments: 
##-- population map



#########################


## Generate list of populations, and number of samples per population
import sys
popmap = open(sys.argv[1], "r")
script = open("Eleni_filter_by_MinorAlleleFrequency_takeARGS.py", "w")
replace = open("replace.txt", "r")
replace_str = replace.read()
replace.close()

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



## Start writing the MAF filtering script ##

script.write("### This scripts corrects for minor allele frequency\n#Adjusted from Eleni's script (June 15,2015) to take arguments\n#MF 2/24/2017\n\n#################\n\n\nimport sys\n\n")

script.write("# Open your files for reading and writing\ngenotypes_file = open(sys.argv[1],'r')\noutput_freqs = open(sys.argv[2],'w')\nfiltered_genotypes = open(sys.argv[3],'w')\nblacklisted_genotypes = open(sys.argv[4],'w')\nblacklisted_MAF = open(sys.argv[5],'w')\n\n# Tell the computer that your files have headers\nheader = True\n\n")


# adding headers to files

script.write("# This code creates a list of each allele for each population. This will be the headers for the file that outputs the allele frequencies. Modify as needed.\n")

MAFhead = "Locus"
for pop in PopList: 
	newstr = "\t" + "Allele1_" + pop + "\t" + "Allele2_" + pop
	MAFhead = MAFhead + newstr

script.write("MAF_header = '" + MAFhead + "'\n")
script.write(r"output_freqs.write(MAF_header + '\n')") #"r" makes it a raw string
script.write("\n" + r"blacklisted_MAF.write(MAF_header + '\n')")
script.write("\n\n\n\n")


# create column counts
column = 1
column_indices = {}
for pop in PopList: 
	n_samples = PopDict[pop]
	new_column = column + n_samples
	column_indices[pop] = [column, new_column]
	column = new_column

	



# start "for mystring in genotypes_file" loop

script.write("for mystring in genotypes_file:\n\tif header:\n\t\tgenotypes_header = mystring\n\t\tfiltered_genotypes.write(genotypes_header)\n\t\tblacklisted_genotypes.write(genotypes_header)\n\t\theader = False\n\telse:\n\t\t" + r"stripped_string = mystring.strip('\n')" + "\n\t\tlocus = stripped_string.split(',')[0]\n\t\tlocus_freqs = []\n\t\tbad_locus_freqs = []\n") 

# generate a string that splits by columns for each population and write to script

for pop in PopList: 
	cols = column_indices[pop]
	newstr = pop + " = stripped_string.split(',')[" + str(cols[0]) + ":" + str(cols[1]) + "]"
	script.write("\t\t" + newstr + "\n")
script.write("\n\n")


# generate a string that counts homozygotes and heterozygotes

for pop in PopList: 
	newstr = "\t\tCountOf_homo1_" + pop + " = float(" + pop + ".count('0101'))\n"
	newstr = newstr + "\t\tCountOf_homo2_" + pop + " = float(" + pop + ".count('0202'))\n"
	newstr = newstr + "\t\tCountOf_het_" + pop + " = float(" + pop + ".count('0102'))\n"
	script.write(newstr + "\n\n")

# generate a string that calculates allele frequencies for each population. 

for pop in PopList: 
	newstr = "\t\t" + "total_alleles_" + pop + "=2*(CountOf_homo1_" + pop + " + CountOf_homo2_" + pop + " + CountOf_het_" + pop + " + 0.000000001)" + "\n\t\t" + "FrequencyOf_allele1_" + pop + " = ((2 * CountOf_homo1_" + pop + ") + (CountOf_het_" + pop + ")) / (total_alleles_" + pop + ")" + "\n\t\t" + "FrequencyOf_allele2_" + pop + " = ((2 * CountOf_homo2_" + pop + ") + (CountOf_het_" + pop + ")) / (total_alleles_" + pop + ")\n\n"
	script.write(newstr)


script.write("\t\tif (")
reps = 1
for pop in PopList:
	if reps == 1: 
		newstr = "(FrequencyOf_allele1_" + pop + " >= 0.05)"
		script.write(newstr)
	else: 
		newstr = " or (FrequencyOf_allele1_" + pop + " >= 0.05)"
		script.write(newstr)
	reps += 1

script.write(") and (")
reps = 1
for pop in PopList: 
	if reps == 1: 
		newstr = "(FrequencyOf_allele2_" + pop + " >= 0.05)"
		script.write(newstr)
	else: 
		newstr = " or (FrequencyOf_allele2_" + pop + " >= 0.05)"
		script.write(newstr)
	reps += 1
script.write("):\n")

script.write("\t\t\t" + r"locus_freqs.append(locus+'\t'+")
end = len(PopList)
reps2 = 1
for pop in PopList: 
	if reps2 >= end:
		newstr = "str(FrequencyOf_allele1_" + pop + ")"
		script.write(newstr)
	else:
		newstr = "str(FrequencyOf_allele1_" + pop + r") + '\t' + "
		script.write(newstr)
	reps2 += 1
script.write(")\n\n")

script.write("\t\t\t" + r"locus_write = str(locus_freqs)." + replace_str + "\n")

script.write("\t\t\t" + r"output_freqs.write(locus_write + '\n')" + "\n\t\t\tfiltered_genotypes.write(mystring)\n")
script.write("\t\telse:\n")

script.write("\t\t\t" + r"bad_locus_freqs.append(locus+'\t'+ ")
end = len(PopList)
reps3 = 1
for pop in PopList: 
	if reps3 >= end:
		newstr = "str(FrequencyOf_allele1_" + pop+ ")"
	else:
		newstr = "str(FrequencyOf_allele1_" + pop+ r") + '\t' + "
	script.write(newstr)
	reps3 += 1
script.write(")\n")

script.write("\t\t\t" + r"bad_locus_write = str(bad_locus_freqs)." + replace_str + "\n")
script.write("\t\t\tprint bad_locus_write\n" + "\t\t\t" + r"blacklisted_MAF.write(bad_locus_write + '\n')" + "\n\t\t\tblacklisted_genotypes.write(mystring)\n\n")

script.write("#close open files\ngenotypes_file.close()\nblacklisted_genotypes.close()\nblacklisted_MAF.close()\nfiltered_genotypes.close()\noutput_freqs.close()")
popmap.close()
script.close()
