### This scripts corrects for minor allele frequency
#Adjusted from Eleni's script (June 15,2015) to take arguments
#MF 2/24/2017

#################


import sys

# Open your files for reading and writing
genotypes_file = open(sys.argv[1],'r')
output_freqs = open(sys.argv[2],'w')
filtered_genotypes = open(sys.argv[3],'w')
blacklisted_genotypes = open(sys.argv[4],'w')
blacklisted_MAF = open(sys.argv[5],'w')

# Tell the computer that your files have headers
header = True

# This code creates a list of each allele for each population. This will be the headers for the file that outputs the allele frequencies. Modify as needed.
MAF_header = 'Locus	Allele1_Pohang15	Allele2_Pohang15	Allele1_Geoje15	Allele2_Geoje15	Allele1_Namhae15	Allele2_Namhae15	Allele1_YellowSea16	Allele2_YellowSea16	Allele1_Jukbyeon07	Allele2_Jukbyeon07	Allele1_JinhaeBay07	Allele2_JinhaeBay07	Allele1_JinhaeBay08	Allele2_JinhaeBay08	Allele1_Boryeong07	Allele2_Boryeong07	Allele1_Geoje14	Allele2_Geoje14'
output_freqs.write(MAF_header + '\n')
blacklisted_MAF.write(MAF_header + '\n')



for mystring in genotypes_file:
	if header:
		genotypes_header = mystring
		filtered_genotypes.write(genotypes_header)
		blacklisted_genotypes.write(genotypes_header)
		header = False
	else:
		stripped_string = mystring.strip('\n')
		locus = stripped_string.split(',')[0]
		locus_freqs = []
		bad_locus_freqs = []
		Pohang15 = stripped_string.split(',')[1:32]
		Geoje15 = stripped_string.split(',')[32:65]
		Namhae15 = stripped_string.split(',')[65:81]
		YellowSea16 = stripped_string.split(',')[81:105]
		Jukbyeon07 = stripped_string.split(',')[105:121] + stripped_string.split(',')[122:130]
		JinhaeBay07 = stripped_string.split(',')[130:132] + stripped_string.split(',')[134:134] + stripped_string.split(',')[136:138] + stripped_string.split(',')[139:139] + stripped_string.split(',')[141:141] + stripped_string.split(',')[143:146] + stripped_string.split(',')[147:147] + stripped_string.split(',')[150:150] + stripped_string.split(',')[153:157] + stripped_string.split(',')[158:160] + stripped_string.split(',')[161:163] + stripped_string.split(',')[164:164]
		JinhaeBay08 = stripped_string.split(',')[167:167] + stripped_string.split(',')[169:175] + stripped_string.split(',')[177:179] + stripped_string.split(',')[180:185] + stripped_string.split(',')[186:186] + stripped_string.split(',')[188:190] + stripped_string.split(',')[191:191] + stripped_string.split(',')[193:193] + stripped_string.split(',')[195:200] + stripped_string.split(',')[201:201] + stripped_string.split(',')[203:210]
		Boryeong07 = stripped_string.split(',')[211:230] + stripped_string.split(',')[231:231]
		Geoje14 = stripped_string.split(',')[233:235] + stripped_string.split(',')[236:256] + stripped_string.split(',')[257:265]


		CountOf_homo1_Pohang15 = float(Pohang15.count('0101'))
		CountOf_homo2_Pohang15 = float(Pohang15.count('0202'))
		CountOf_het_Pohang15 = float(Pohang15.count('0102'))


		CountOf_homo1_Geoje15 = float(Geoje15.count('0101'))
		CountOf_homo2_Geoje15 = float(Geoje15.count('0202'))
		CountOf_het_Geoje15 = float(Geoje15.count('0102'))


		CountOf_homo1_Namhae15 = float(Namhae15.count('0101'))
		CountOf_homo2_Namhae15 = float(Namhae15.count('0202'))
		CountOf_het_Namhae15 = float(Namhae15.count('0102'))


		CountOf_homo1_YellowSea16 = float(YellowSea16.count('0101'))
		CountOf_homo2_YellowSea16 = float(YellowSea16.count('0202'))
		CountOf_het_YellowSea16 = float(YellowSea16.count('0102'))


		CountOf_homo1_Jukbyeon07 = float(Jukbyeon07.count('0101'))
		CountOf_homo2_Jukbyeon07 = float(Jukbyeon07.count('0202'))
		CountOf_het_Jukbyeon07 = float(Jukbyeon07.count('0102'))


		CountOf_homo1_JinhaeBay07 = float(JinhaeBay07.count('0101'))
		CountOf_homo2_JinhaeBay07 = float(JinhaeBay07.count('0202'))
		CountOf_het_JinhaeBay07 = float(JinhaeBay07.count('0102'))


		CountOf_homo1_JinhaeBay08 = float(JinhaeBay08.count('0101'))
		CountOf_homo2_JinhaeBay08 = float(JinhaeBay08.count('0202'))
		CountOf_het_JinhaeBay08 = float(JinhaeBay08.count('0102'))


		CountOf_homo1_Boryeong07 = float(Boryeong07.count('0101'))
		CountOf_homo2_Boryeong07 = float(Boryeong07.count('0202'))
		CountOf_het_Boryeong07 = float(Boryeong07.count('0102'))


		CountOf_homo1_Geoje14 = float(Geoje14.count('0101'))
		CountOf_homo2_Geoje14 = float(Geoje14.count('0202'))
		CountOf_het_Geoje14 = float(Geoje14.count('0102'))


		total_alleles_Pohang15=2*(CountOf_homo1_Pohang15 + CountOf_homo2_Pohang15 + CountOf_het_Pohang15 + 0.000000001)
		FrequencyOf_allele1_Pohang15 = ((2 * CountOf_homo1_Pohang15) + (CountOf_het_Pohang15)) / (total_alleles_Pohang15)
		FrequencyOf_allele2_Pohang15 = ((2 * CountOf_homo2_Pohang15) + (CountOf_het_Pohang15)) / (total_alleles_Pohang15)

		total_alleles_Geoje15=2*(CountOf_homo1_Geoje15 + CountOf_homo2_Geoje15 + CountOf_het_Geoje15 + 0.000000001)
		FrequencyOf_allele1_Geoje15 = ((2 * CountOf_homo1_Geoje15) + (CountOf_het_Geoje15)) / (total_alleles_Geoje15)
		FrequencyOf_allele2_Geoje15 = ((2 * CountOf_homo2_Geoje15) + (CountOf_het_Geoje15)) / (total_alleles_Geoje15)

		total_alleles_Namhae15=2*(CountOf_homo1_Namhae15 + CountOf_homo2_Namhae15 + CountOf_het_Namhae15 + 0.000000001)
		FrequencyOf_allele1_Namhae15 = ((2 * CountOf_homo1_Namhae15) + (CountOf_het_Namhae15)) / (total_alleles_Namhae15)
		FrequencyOf_allele2_Namhae15 = ((2 * CountOf_homo2_Namhae15) + (CountOf_het_Namhae15)) / (total_alleles_Namhae15)

		total_alleles_YellowSea16=2*(CountOf_homo1_YellowSea16 + CountOf_homo2_YellowSea16 + CountOf_het_YellowSea16 + 0.000000001)
		FrequencyOf_allele1_YellowSea16 = ((2 * CountOf_homo1_YellowSea16) + (CountOf_het_YellowSea16)) / (total_alleles_YellowSea16)
		FrequencyOf_allele2_YellowSea16 = ((2 * CountOf_homo2_YellowSea16) + (CountOf_het_YellowSea16)) / (total_alleles_YellowSea16)

		total_alleles_Jukbyeon07=2*(CountOf_homo1_Jukbyeon07 + CountOf_homo2_Jukbyeon07 + CountOf_het_Jukbyeon07 + 0.000000001)
		FrequencyOf_allele1_Jukbyeon07 = ((2 * CountOf_homo1_Jukbyeon07) + (CountOf_het_Jukbyeon07)) / (total_alleles_Jukbyeon07)
		FrequencyOf_allele2_Jukbyeon07 = ((2 * CountOf_homo2_Jukbyeon07) + (CountOf_het_Jukbyeon07)) / (total_alleles_Jukbyeon07)

		total_alleles_JinhaeBay07=2*(CountOf_homo1_JinhaeBay07 + CountOf_homo2_JinhaeBay07 + CountOf_het_JinhaeBay07 + 0.000000001)
		FrequencyOf_allele1_JinhaeBay07 = ((2 * CountOf_homo1_JinhaeBay07) + (CountOf_het_JinhaeBay07)) / (total_alleles_JinhaeBay07)
		FrequencyOf_allele2_JinhaeBay07 = ((2 * CountOf_homo2_JinhaeBay07) + (CountOf_het_JinhaeBay07)) / (total_alleles_JinhaeBay07)

		total_alleles_JinhaeBay08=2*(CountOf_homo1_JinhaeBay08 + CountOf_homo2_JinhaeBay08 + CountOf_het_JinhaeBay08 + 0.000000001)
		FrequencyOf_allele1_JinhaeBay08 = ((2 * CountOf_homo1_JinhaeBay08) + (CountOf_het_JinhaeBay08)) / (total_alleles_JinhaeBay08)
		FrequencyOf_allele2_JinhaeBay08 = ((2 * CountOf_homo2_JinhaeBay08) + (CountOf_het_JinhaeBay08)) / (total_alleles_JinhaeBay08)

		total_alleles_Boryeong07=2*(CountOf_homo1_Boryeong07 + CountOf_homo2_Boryeong07 + CountOf_het_Boryeong07 + 0.000000001)
		FrequencyOf_allele1_Boryeong07 = ((2 * CountOf_homo1_Boryeong07) + (CountOf_het_Boryeong07)) / (total_alleles_Boryeong07)
		FrequencyOf_allele2_Boryeong07 = ((2 * CountOf_homo2_Boryeong07) + (CountOf_het_Boryeong07)) / (total_alleles_Boryeong07)

		total_alleles_Geoje14=2*(CountOf_homo1_Geoje14 + CountOf_homo2_Geoje14 + CountOf_het_Geoje14 + 0.000000001)
		FrequencyOf_allele1_Geoje14 = ((2 * CountOf_homo1_Geoje14) + (CountOf_het_Geoje14)) / (total_alleles_Geoje14)
		FrequencyOf_allele2_Geoje14 = ((2 * CountOf_homo2_Geoje14) + (CountOf_het_Geoje14)) / (total_alleles_Geoje14)

		if ((FrequencyOf_allele1_Pohang15 >= 0.05) or (FrequencyOf_allele1_Geoje15 >= 0.05) or (FrequencyOf_allele1_Namhae15 >= 0.05) or (FrequencyOf_allele1_YellowSea16 >= 0.05) or (FrequencyOf_allele1_Jukbyeon07 >= 0.05) or (FrequencyOf_allele1_JinhaeBay07 >= 0.05) or (FrequencyOf_allele1_JinhaeBay08 >= 0.05) or (FrequencyOf_allele1_Boryeong07 >= 0.05) or (FrequencyOf_allele1_Geoje14 >= 0.05)) and ((FrequencyOf_allele2_Pohang15 >= 0.05) or (FrequencyOf_allele2_Geoje15 >= 0.05) or (FrequencyOf_allele2_Namhae15 >= 0.05) or (FrequencyOf_allele2_YellowSea16 >= 0.05) or (FrequencyOf_allele2_Jukbyeon07 >= 0.05) or (FrequencyOf_allele2_JinhaeBay07 >= 0.05) or (FrequencyOf_allele2_JinhaeBay08 >= 0.05) or (FrequencyOf_allele2_Boryeong07 >= 0.05) or (FrequencyOf_allele2_Geoje14 >= 0.05)):
			locus_freqs.append(locus+'\t'+str(FrequencyOf_allele1_Pohang15) + '\t' + str(FrequencyOf_allele1_Geoje15) + '\t' + str(FrequencyOf_allele1_Namhae15) + '\t' + str(FrequencyOf_allele1_YellowSea16) + '\t' + str(FrequencyOf_allele1_Jukbyeon07) + '\t' + str(FrequencyOf_allele1_JinhaeBay07) + '\t' + str(FrequencyOf_allele1_JinhaeBay08) + '\t' + str(FrequencyOf_allele1_Boryeong07) + '\t' + str(FrequencyOf_allele1_Geoje14))

			locus_write = str(locus_freqs).replace('[','').replace(',','\t').replace(']', '').replace("'", '').replace(' ','').replace('\\n','').replace('\\t','\t')

			output_freqs.write(locus_write + '\n')
			filtered_genotypes.write(mystring)
		else:
			bad_locus_freqs.append(locus+'\t'+ str(FrequencyOf_allele1_Pohang15) + '\t' + str(FrequencyOf_allele1_Geoje15) + '\t' + str(FrequencyOf_allele1_Namhae15) + '\t' + str(FrequencyOf_allele1_YellowSea16) + '\t' + str(FrequencyOf_allele1_Jukbyeon07) + '\t' + str(FrequencyOf_allele1_JinhaeBay07) + '\t' + str(FrequencyOf_allele1_JinhaeBay08) + '\t' + str(FrequencyOf_allele1_Boryeong07) + '\t' + str(FrequencyOf_allele1_Geoje14))
			bad_locus_write = str(bad_locus_freqs).replace('[','').replace(',','\t').replace(']', '').replace("'", '').replace(' ','').replace('\\n','').replace('\\t','\t')
			blacklisted_MAF.write(bad_locus_write + '\n')
			blacklisted_genotypes.write(mystring)

#close open files
genotypes_file.close()
blacklisted_genotypes.close()
blacklisted_MAF.close()
filtered_genotypes.close()
output_freqs.close()