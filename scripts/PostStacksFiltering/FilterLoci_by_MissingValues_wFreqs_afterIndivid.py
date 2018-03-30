### This scripts removes loci with too much missing data (you set the threshold)
#Adjusted from Eleni's script (June 15,2015) to take arguments
#MF 2/28/2017

#################


import sys

# Open your files for reading and writing
genotypes_file = open(sys.argv[1],'r')
clean_output_file = open(sys.argv[2],'w')
good_freqs = open(sys.argv[3], 'w')
blacklisted_output_file = open(sys.argv[4], 'w')
bad_freqs = open(sys.argv[5], 'w')
freq = float(sys.argv[6])


#run for loop to counting missing genotypes by locus for each population

count = 0
bad_count = 0

for mystring in genotypes_file:		# Read in each line in the file as a string
	if count == 0:
		genotypes_header = mystring
		clean_output_file.write(genotypes_header)
		blacklisted_output_file.write(genotypes_header)
		pops = ["Pohang15", "Geoje15", "Namhae15", "YellowSea16", "Jukbyeon07", "JinhaeBay07", "JinhaeBay08","Boryeong07","Geoje14","Kodiak03","Adak06","WashCoast05","HecStrait04", "SalishSea12_13","JuandeFuca12","PWSound12","UnimakPass03"]
		good_freqs.write("\t".join(pops)+"\n")
		bad_freqs.write("\t".join(pops)+"\n")
		count += 1
	else:
		count += 1
		overall_percent_missingdata = []
		stripped_string = mystring.strip('\n')
		locus_name = stripped_string.split(',')[0]
		Pohang15 = stripped_string.split(',')[1:35]
		Geoje15 = stripped_string.split(',')[35:72]
		Namhae15 = stripped_string.split(',')[72:91]
		YellowSea16 = stripped_string.split(',')[91:121]
		Jukbyeon07 = stripped_string.split(',')[121:158]
		JinhaeBay07 = stripped_string.split(',')[158:211]
		JinhaeBay08 = stripped_string.split(',')[211:267]
		Boryeong07 = stripped_string.split(',')[267:291]
		Geoje14 = stripped_string.split(',')[291:327]
		Kodiak03 = stripped_string.split(',')[327:375]
		Adak06 = stripped_string.split(',')[375:423]
		WashCoast05 = stripped_string.split(',')[423:471]
		HecStrait04 = stripped_string.split(',')[471:519]
		SalishSea12_13 = stripped_string.split(',')[519:544]
		SalishSea12_13.append(stripped_string.split(',')[567:584])
		JuandeFuca12 = stripped_string.split(',')[544:567]
		PWSound12 = stripped_string.split(',')[584:632]
		UnimakPass03 = stripped_string.split(',')[632:680]

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
#next pop
		Count_MissingGenotypesByLocus_Kodiak03 = float(Kodiak03.count('0000'))
		NumberOf_Kodiak03_individuals = float(len(Kodiak03))
		Percent_MissingData_Kodiak03 = float(Count_MissingGenotypesByLocus_Kodiak03/NumberOf_Kodiak03_individuals)
		overall_percent_missingdata.append(Percent_MissingData_Kodiak03)
#next pop
		Count_MissingGenotypesByLocus_Adak06 = float(Adak06.count('0000'))
		NumberOf_Adak06_individuals = float(len(Adak06))
		Percent_MissingData_Adak06 = float(Count_MissingGenotypesByLocus_Adak06/NumberOf_Adak06_individuals)
		overall_percent_missingdata.append(Percent_MissingData_Adak06)
#next pop
		Count_MissingGenotypesByLocus_WashCoast05 = float(WashCoast05.count('0000'))
		NumberOf_WashCoast05_individuals = float(len(WashCoast05))
		Percent_MissingData_WashCoast05 = float(Count_MissingGenotypesByLocus_WashCoast05/NumberOf_WashCoast05_individuals)
		overall_percent_missingdata.append(Percent_MissingData_WashCoast05)
#next pop
		Count_MissingGenotypesByLocus_HecStrait04 = float(HecStrait04.count('0000'))
		NumberOf_HecStrait04_individuals = float(len(HecStrait04))
		Percent_MissingData_HecStrait04 = float(Count_MissingGenotypesByLocus_HecStrait04/NumberOf_HecStrait04_individuals)
		overall_percent_missingdata.append(Percent_MissingData_HecStrait04)
#next pop
		Count_MissingGenotypesByLocus_SalishSea12_13 = float(SalishSea12_13.count('0000'))
		NumberOf_SalishSea12_13_individuals = float(len(SalishSea12_13))
		Percent_MissingData_SalishSea12_13 = float(Count_MissingGenotypesByLocus_SalishSea12_13/NumberOf_SalishSea12_13_individuals)
		overall_percent_missingdata.append(Percent_MissingData_SalishSea12_13)
#next pop
		Count_MissingGenotypesByLocus_JuandeFuca12 = float(JuandeFuca12.count('0000'))
		NumberOf_JuandeFuca12_individuals = float(len(JuandeFuca12))
		Percent_MissingData_JuandeFuca12 = float(Count_MissingGenotypesByLocus_JuandeFuca12/NumberOf_JuandeFuca12_individuals)
		overall_percent_missingdata.append(Percent_MissingData_JuandeFuca12)
#next pop
		Count_MissingGenotypesByLocus_PWSound12 = float(PWSound12.count('0000'))
		NumberOf_PWSound12_individuals = float(len(PWSound12))
		Percent_MissingData_PWSound12 = float(Count_MissingGenotypesByLocus_PWSound12/NumberOf_PWSound12_individuals)
		overall_percent_missingdata.append(Percent_MissingData_PWSound12)
#next pop
		Count_MissingGenotypesByLocus_UnimakPass03 = float(UnimakPass03.count('0000'))
		NumberOf_UnimakPass03_individuals = float(len(UnimakPass03))
		Percent_MissingData_UnimakPass03 = float(Count_MissingGenotypesByLocus_UnimakPass03/NumberOf_UnimakPass03_individuals)
		overall_percent_missingdata.append(Percent_MissingData_UnimakPass03)

#write loci to appropriate file
		if all(i < freq for i in overall_percent_missingdata):
			clean_output_file.write(mystring)
			good_freqs.write(locus_name + "\t" + str(overall_percent_missingdata) + "\n")
		else: 
			blacklisted_output_file.write(mystring)
			bad_freqs.write(locus_name + "\t" + str(overall_percent_missingdata) + "\n")
			bad_count += 1
#print output
n_loci = str(count - 1)
print 'processed ' + n_loci + ' loci'
print 'Number of loci removed: ' + str(bad_count)