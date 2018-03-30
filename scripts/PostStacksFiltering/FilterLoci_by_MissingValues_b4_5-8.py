### This scripts removes loci with too much missing data (you set the threshold)
#Adjusted from Eleni's script (June 15,2015) to take arguments
#MF 2/28/2017

#################


import sys

# Open your files for reading and writing
genotypes_file = open(sys.argv[1],'r')
clean_output_file = open(sys.argv[2],'w')
blacklisted_output_file = open(sys.argv[3], 'w')


#run for loop to counting missing genotypes by locus for each population

count = 0
bad_count = 0

for mystring in genotypes_file:		# Read in each line in the file as a string
	if count == 0:
		genotypes_header = mystring
		clean_output_file.write(genotypes_header)
		blacklisted_output_file.write(genotypes_header)
		count += 1
	else:
		count += 1
		overall_percent_missingdata = []
		stripped_string = mystring.strip('\n')
		locus_name = stripped_string.split(',')[0]
		Pohang15 = stripped_string.split(',')[1:32]
		Geoje15 = stripped_string.split(',')[32:65]
		Namhae15 = stripped_string.split(',')[65:81]
		YellowSea16 = stripped_string.split(',')[81:105]
		Jukbyeon07 = stripped_string.split(',')[105:121] + stripped_string.split(',')[122:130]
		JinhaeBay07 = stripped_string.split(',')[130:132] + stripped_string.split(',')[134:134] + stripped_string.split(',')[136:138] + stripped_string.split(',')[139:139] + stripped_string.split(',')[141:141] + stripped_string.split(',')[143:146] + stripped_string.split(',')[147:147] + stripped_string.split(',')[150:150] + stripped_string.split(',')[153:157] + stripped_string.split(',')[158:160] + stripped_string.split(',')[161:163] + stripped_string.split(',')[164:164]
		JinhaeBay08 = stripped_string.split(',')[167:167] + stripped_string.split(',')[169:175] + stripped_string.split(',')[177:179] + stripped_string.split(',')[180:185] + stripped_string.split(',')[186:186] + stripped_string.split(',')[188:190] + stripped_string.split(',')[191:191] + stripped_string.split(',')[193:193] + stripped_string.split(',')[195:200] + stripped_string.split(',')[201:201] + stripped_string.split(',')[203:210]
		Boryeong07 = stripped_string.split(',')[211:230] + stripped_string.split(',')[231:231]
		Geoje14 = stripped_string.split(',')[233:235] + stripped_string.split(',')[236:256] + stripped_string.split(',')[257:265]


#per population counts
#next pop
		Count_MissingGenotypesByLocus_Pohang15 = float(Pohang15.count('0000'))
		NumberOf_Pohang15_individuals = float(len(Pohang15))
		Percent_MissingData_Pohang15 = float(Count_MissingGenotypesByLocus_Pohang15/NumberOf_Pohang15_individuals)
		overall_percent_missingdata.append(Percent_MissingData_Pohang15)
#next pop
		Count_MissingGenotypesByLocus_Geoje15 = float(Geoje15.count('0000'))
		NumberOf_Geoje15_individuals = float(len(Geoje15))
		Percent_MissingData_Geoje15 = float(Count_MissingGenotypesByLocus_Geoje15/NumberOf_Geoje15_individuals)
		overall_percent_missingdata.append(Percent_MissingData_Geoje15)
#next pop
		Count_MissingGenotypesByLocus_Namhae15 = float(Namhae15.count('0000'))
		NumberOf_Namhae15_individuals = float(len(Namhae15))
		Percent_MissingData_Namhae15 = float(Count_MissingGenotypesByLocus_Namhae15/NumberOf_Namhae15_individuals)
		overall_percent_missingdata.append(Percent_MissingData_Namhae15)
#next pop
		Count_MissingGenotypesByLocus_YellowSea16 = float(YellowSea16.count('0000'))
		NumberOf_YellowSea16_individuals = float(len(YellowSea16))
		Percent_MissingData_YellowSea16 = float(Count_MissingGenotypesByLocus_YellowSea16/NumberOf_YellowSea16_individuals)
		overall_percent_missingdata.append(Percent_MissingData_YellowSea16)
#next pop
		Count_MissingGenotypesByLocus_Jukbyeon07 = float(Jukbyeon07.count('0000'))
		NumberOf_Jukbyeon07_individuals = float(len(Jukbyeon07))
		Percent_MissingData_Jukbyeon07 = float(Count_MissingGenotypesByLocus_Jukbyeon07/NumberOf_Jukbyeon07_individuals)
		overall_percent_missingdata.append(Percent_MissingData_Jukbyeon07)
#next pop
		Count_MissingGenotypesByLocus_JinhaeBay07 = float(JinhaeBay07.count('0000'))
		NumberOf_JinhaeBay07_individuals = float(len(JinhaeBay07))
		Percent_MissingData_JinhaeBay07 = float(Count_MissingGenotypesByLocus_JinhaeBay07/NumberOf_JinhaeBay07_individuals)
		overall_percent_missingdata.append(Percent_MissingData_JinhaeBay07)
#next pop
		Count_MissingGenotypesByLocus_JinhaeBay08 = float(JinhaeBay08.count('0000'))
		NumberOf_JinhaeBay08_individuals = float(len(JinhaeBay08))
		Percent_MissingData_JinhaeBay08 = float(Count_MissingGenotypesByLocus_JinhaeBay08/NumberOf_JinhaeBay08_individuals)
		overall_percent_missingdata.append(Percent_MissingData_JinhaeBay08)
#next pop
		Count_MissingGenotypesByLocus_Boryeong07 = float(Boryeong07.count('0000'))
		NumberOf_Boryeong07_individuals = float(len(Boryeong07))
		Percent_MissingData_Boryeong07 = float(Count_MissingGenotypesByLocus_Boryeong07/NumberOf_Boryeong07_individuals)
		overall_percent_missingdata.append(Percent_MissingData_Boryeong07)
#next pop
		Count_MissingGenotypesByLocus_Geoje14 = float(Geoje14.count('0000'))
		NumberOf_Geoje14_individuals = float(len(Geoje14))
		Percent_MissingData_Geoje14 = float(Count_MissingGenotypesByLocus_Geoje14/NumberOf_Geoje14_individuals)
		overall_percent_missingdata.append(Percent_MissingData_Geoje14)

#write loci to appropriate file
		if all(i < 0.50 for i in overall_percent_missingdata):
			clean_output_file.write(mystring)
		else: 
			blacklisted_output_file.write(mystring)
			bad_count += 1
#print output
n_loci = str(count - 1)
print 'processed ' + n_loci + ' loci'
print 'Number of loci removed: ' + str(bad_count)