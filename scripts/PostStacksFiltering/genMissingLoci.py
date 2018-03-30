### This file generates Eleni's Missing Loci filtering file ###

## MF 2/24/2017

## Arguments: 
##-- population map



#########################

import sys
popmap = open(sys.argv[1], "r")
script = open("FilterLoci_by_MissingValues.py", "w")


## Start writing the MAF filtering script ##

script.write("### This scripts removes loci with too much missing data (you set the threshold)\n#Adjusted from Eleni's script (June 15,2015) to take arguments\n#MF 2/28/2017\n\n#################\n\n\nimport sys\n\n")


script.write("# Open your files for reading and writing\ngenotypes_file = open(sys.argv[1],'r')\nclean_output_file = open(sys.argv[2],'w')\nblacklisted_output_file = open(sys.argv[3], 'w')\n\n")


## Count the missing genotypes in each population

#--- generate the same ordered dictionary used in the MAF filtering
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


#--- initiate the for loop
script.write("\n#run for loop to counting missing genotypes by locus for each population\n\ncount = 0\nbad_count = 0\n\nfor mystring in genotypes_file:		# Read in each line in the file as a string\n")
script.write("\tif count == 0:\n\t\tgenotypes_header = mystring\n\t\tclean_output_file.write(genotypes_header)\n\t\tblacklisted_output_file.write(genotypes_header)\n\t\tcount += 1\n")

script.write("\telse:\n\t\tcount += 1\n\t\toverall_percent_missingdata = []\n\t\t" + r"stripped_string = mystring.strip('\n')" + "\n\t\t" + "locus_name = stripped_string.split(',')[0]" + "\n")

#--- use dictionary to generate column indices (should be the same order as in MAF filtering)
PopList = PopDict.keys()
column = 1
column_indices = {}
for pop in PopList: 
	n_samples = PopDict[pop]
	new_column = column + n_samples
	column_indices[pop] = [column, new_column]
	column = new_column


for pop in PopList: 
	cols = column_indices[pop]
	newstr = pop + " = stripped_string.split(',')[" + str(cols[0]) + ":" + str(cols[1]) + "]"
	script.write("\t\t" + newstr + "\n")
script.write("\n#per population counts\n")


#--- count missing data in each population
for pop in PopList: 
	newstr = "#next pop" + "\n\t\tCount_MissingGenotypesByLocus_" + pop + " = float(" + pop + ".count('0000'))" + "\n\t\tNumberOf_" + pop + "_individuals = float(len(" + pop + "))\n\t\tPercent_MissingData_" + pop + " = float(Count_MissingGenotypesByLocus_" + pop + "/NumberOf_" + pop + "_individuals)\n\t\t" + "overall_percent_missingdata.append(Percent_MissingData_" + pop + ")\n"
	script.write(newstr)

#--- write good loci to one file, bad loci to another

script.write("\n#write loci to appropriate file\n\t\tif all(i < 0.50 for i in overall_percent_missingdata):\n\t\t\tclean_output_file.write(mystring)\n\t\telse: \n\t\t\tblacklisted_output_file.write(mystring)\n\t\t\tbad_count += 1")

#--- print output when finished

script.write("\n#print output\nn_loci = str(count - 1)\n" + "print 'processed ' + n_loci + ' loci'" + "\n" + "print 'Number of loci removed: ' + str(bad_count)")

script.close()
popmap.close()
