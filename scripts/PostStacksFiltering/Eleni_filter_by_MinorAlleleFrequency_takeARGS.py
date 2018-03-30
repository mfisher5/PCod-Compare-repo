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
MAF_header = 'Locus	Allele1_Pohang15	Allele2_Pohang15	Allele1_Geoje15	Allele2_Geoje15	Allele1_Namhae15	Allele2_Namhae15	Allele1_YellowSea16	Allele2_YellowSea16	Allele1_Jukbyeon07	Allele2_Jukbyeon07	Allele1_JinhaeBay07	Allele2_JinhaeBay07	Allele1_JinhaeBay08	Allele2_JinhaeBay08	Allele1_Boryeong07	Allele2_Boryeong07	Allele1_Geoje14	Allele2_Geoje14	Allele1_Kodiak03	Allele2_Kodiak03	Allele1_Adak06	Allele2_Adak06	Allele1_WashCoast05	Allele2_WashCoast05	Allele1_HecStrait04	Allele2_HecStrait04	Allele1_SalishSea12_13	Allele2_SalishSea12_13	Allele1_JuandeFuca12	Allele2_JuandeFuca12	Allele1_PWSound12	Allele2_PWSound12	Allele1_UnimakPass03	Allele2_UnimakPass03'
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
		SalishSea12_13 = stripped_string.split(',')[519:561]
		JuandeFuca12 = stripped_string.split(',')[561:584]
		PWSound12 = stripped_string.split(',')[584:632]
		UnimakPass03 = stripped_string.split(',')[632:680]


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


		CountOf_homo1_Kodiak03 = float(Kodiak03.count('0101'))
		CountOf_homo2_Kodiak03 = float(Kodiak03.count('0202'))
		CountOf_het_Kodiak03 = float(Kodiak03.count('0102'))


		CountOf_homo1_Adak06 = float(Adak06.count('0101'))
		CountOf_homo2_Adak06 = float(Adak06.count('0202'))
		CountOf_het_Adak06 = float(Adak06.count('0102'))


		CountOf_homo1_WashCoast05 = float(WashCoast05.count('0101'))
		CountOf_homo2_WashCoast05 = float(WashCoast05.count('0202'))
		CountOf_het_WashCoast05 = float(WashCoast05.count('0102'))


		CountOf_homo1_HecStrait04 = float(HecStrait04.count('0101'))
		CountOf_homo2_HecStrait04 = float(HecStrait04.count('0202'))
		CountOf_het_HecStrait04 = float(HecStrait04.count('0102'))


		CountOf_homo1_SalishSea12_13 = float(SalishSea12_13.count('0101'))
		CountOf_homo2_SalishSea12_13 = float(SalishSea12_13.count('0202'))
		CountOf_het_SalishSea12_13 = float(SalishSea12_13.count('0102'))


		CountOf_homo1_JuandeFuca12 = float(JuandeFuca12.count('0101'))
		CountOf_homo2_JuandeFuca12 = float(JuandeFuca12.count('0202'))
		CountOf_het_JuandeFuca12 = float(JuandeFuca12.count('0102'))


		CountOf_homo1_PWSound12 = float(PWSound12.count('0101'))
		CountOf_homo2_PWSound12 = float(PWSound12.count('0202'))
		CountOf_het_PWSound12 = float(PWSound12.count('0102'))


		CountOf_homo1_UnimakPass03 = float(UnimakPass03.count('0101'))
		CountOf_homo2_UnimakPass03 = float(UnimakPass03.count('0202'))
		CountOf_het_UnimakPass03 = float(UnimakPass03.count('0102'))


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

		total_alleles_Kodiak03=2*(CountOf_homo1_Kodiak03 + CountOf_homo2_Kodiak03 + CountOf_het_Kodiak03 + 0.000000001)
		FrequencyOf_allele1_Kodiak03 = ((2 * CountOf_homo1_Kodiak03) + (CountOf_het_Kodiak03)) / (total_alleles_Kodiak03)
		FrequencyOf_allele2_Kodiak03 = ((2 * CountOf_homo2_Kodiak03) + (CountOf_het_Kodiak03)) / (total_alleles_Kodiak03)

		total_alleles_Adak06=2*(CountOf_homo1_Adak06 + CountOf_homo2_Adak06 + CountOf_het_Adak06 + 0.000000001)
		FrequencyOf_allele1_Adak06 = ((2 * CountOf_homo1_Adak06) + (CountOf_het_Adak06)) / (total_alleles_Adak06)
		FrequencyOf_allele2_Adak06 = ((2 * CountOf_homo2_Adak06) + (CountOf_het_Adak06)) / (total_alleles_Adak06)

		total_alleles_WashCoast05=2*(CountOf_homo1_WashCoast05 + CountOf_homo2_WashCoast05 + CountOf_het_WashCoast05 + 0.000000001)
		FrequencyOf_allele1_WashCoast05 = ((2 * CountOf_homo1_WashCoast05) + (CountOf_het_WashCoast05)) / (total_alleles_WashCoast05)
		FrequencyOf_allele2_WashCoast05 = ((2 * CountOf_homo2_WashCoast05) + (CountOf_het_WashCoast05)) / (total_alleles_WashCoast05)

		total_alleles_HecStrait04=2*(CountOf_homo1_HecStrait04 + CountOf_homo2_HecStrait04 + CountOf_het_HecStrait04 + 0.000000001)
		FrequencyOf_allele1_HecStrait04 = ((2 * CountOf_homo1_HecStrait04) + (CountOf_het_HecStrait04)) / (total_alleles_HecStrait04)
		FrequencyOf_allele2_HecStrait04 = ((2 * CountOf_homo2_HecStrait04) + (CountOf_het_HecStrait04)) / (total_alleles_HecStrait04)

		total_alleles_SalishSea12_13=2*(CountOf_homo1_SalishSea12_13 + CountOf_homo2_SalishSea12_13 + CountOf_het_SalishSea12_13 + 0.000000001)
		FrequencyOf_allele1_SalishSea12_13 = ((2 * CountOf_homo1_SalishSea12_13) + (CountOf_het_SalishSea12_13)) / (total_alleles_SalishSea12_13)
		FrequencyOf_allele2_SalishSea12_13 = ((2 * CountOf_homo2_SalishSea12_13) + (CountOf_het_SalishSea12_13)) / (total_alleles_SalishSea12_13)

		total_alleles_JuandeFuca12=2*(CountOf_homo1_JuandeFuca12 + CountOf_homo2_JuandeFuca12 + CountOf_het_JuandeFuca12 + 0.000000001)
		FrequencyOf_allele1_JuandeFuca12 = ((2 * CountOf_homo1_JuandeFuca12) + (CountOf_het_JuandeFuca12)) / (total_alleles_JuandeFuca12)
		FrequencyOf_allele2_JuandeFuca12 = ((2 * CountOf_homo2_JuandeFuca12) + (CountOf_het_JuandeFuca12)) / (total_alleles_JuandeFuca12)

		total_alleles_PWSound12=2*(CountOf_homo1_PWSound12 + CountOf_homo2_PWSound12 + CountOf_het_PWSound12 + 0.000000001)
		FrequencyOf_allele1_PWSound12 = ((2 * CountOf_homo1_PWSound12) + (CountOf_het_PWSound12)) / (total_alleles_PWSound12)
		FrequencyOf_allele2_PWSound12 = ((2 * CountOf_homo2_PWSound12) + (CountOf_het_PWSound12)) / (total_alleles_PWSound12)

		total_alleles_UnimakPass03=2*(CountOf_homo1_UnimakPass03 + CountOf_homo2_UnimakPass03 + CountOf_het_UnimakPass03 + 0.000000001)
		FrequencyOf_allele1_UnimakPass03 = ((2 * CountOf_homo1_UnimakPass03) + (CountOf_het_UnimakPass03)) / (total_alleles_UnimakPass03)
		FrequencyOf_allele2_UnimakPass03 = ((2 * CountOf_homo2_UnimakPass03) + (CountOf_het_UnimakPass03)) / (total_alleles_UnimakPass03)

		if ((FrequencyOf_allele1_Pohang15 >= 0.05) or (FrequencyOf_allele1_Geoje15 >= 0.05) or (FrequencyOf_allele1_Namhae15 >= 0.05) or (FrequencyOf_allele1_YellowSea16 >= 0.05) or (FrequencyOf_allele1_Jukbyeon07 >= 0.05) or (FrequencyOf_allele1_JinhaeBay07 >= 0.05) or (FrequencyOf_allele1_JinhaeBay08 >= 0.05) or (FrequencyOf_allele1_Boryeong07 >= 0.05) or (FrequencyOf_allele1_Geoje14 >= 0.05) or (FrequencyOf_allele1_Kodiak03 >= 0.05) or (FrequencyOf_allele1_Adak06 >= 0.05) or (FrequencyOf_allele1_WashCoast05 >= 0.05) or (FrequencyOf_allele1_HecStrait04 >= 0.05) or (FrequencyOf_allele1_SalishSea12_13 >= 0.05) or (FrequencyOf_allele1_JuandeFuca12 >= 0.05) or (FrequencyOf_allele1_PWSound12 >= 0.05) or (FrequencyOf_allele1_UnimakPass03 >= 0.05)) and ((FrequencyOf_allele2_Pohang15 >= 0.05) or (FrequencyOf_allele2_Geoje15 >= 0.05) or (FrequencyOf_allele2_Namhae15 >= 0.05) or (FrequencyOf_allele2_YellowSea16 >= 0.05) or (FrequencyOf_allele2_Jukbyeon07 >= 0.05) or (FrequencyOf_allele2_JinhaeBay07 >= 0.05) or (FrequencyOf_allele2_JinhaeBay08 >= 0.05) or (FrequencyOf_allele2_Boryeong07 >= 0.05) or (FrequencyOf_allele2_Geoje14 >= 0.05) or (FrequencyOf_allele2_Kodiak03 >= 0.05) or (FrequencyOf_allele2_Adak06 >= 0.05) or (FrequencyOf_allele2_WashCoast05 >= 0.05) or (FrequencyOf_allele2_HecStrait04 >= 0.05) or (FrequencyOf_allele2_SalishSea12_13 >= 0.05) or (FrequencyOf_allele2_JuandeFuca12 >= 0.05) or (FrequencyOf_allele2_PWSound12 >= 0.05) or (FrequencyOf_allele2_UnimakPass03 >= 0.05)):
			locus_freqs.append(locus+'\t'+str(FrequencyOf_allele1_Pohang15) + '\t' + str(FrequencyOf_allele1_Geoje15) + '\t' + str(FrequencyOf_allele1_Namhae15) + '\t' + str(FrequencyOf_allele1_YellowSea16) + '\t' + str(FrequencyOf_allele1_Jukbyeon07) + '\t' + str(FrequencyOf_allele1_JinhaeBay07) + '\t' + str(FrequencyOf_allele1_JinhaeBay08) + '\t' + str(FrequencyOf_allele1_Boryeong07) + '\t' + str(FrequencyOf_allele1_Geoje14) + '\t' + str(FrequencyOf_allele1_Kodiak03) + '\t' + str(FrequencyOf_allele1_Adak06) + '\t' + str(FrequencyOf_allele1_WashCoast05) + '\t' + str(FrequencyOf_allele1_HecStrait04) + '\t' + str(FrequencyOf_allele1_SalishSea12_13) + '\t' + str(FrequencyOf_allele1_JuandeFuca12) + '\t' + str(FrequencyOf_allele1_PWSound12) + '\t' + str(FrequencyOf_allele1_UnimakPass03))

			locus_write = str(locus_freqs).replace('[','').replace(',','\t').replace(']', '').replace("'", '').replace(' ','').replace('\\n','').replace('\\t','\t')

			output_freqs.write(locus_write + '\n')
			filtered_genotypes.write(mystring)
		else:
			bad_locus_freqs.append(locus+'\t'+ str(FrequencyOf_allele1_Pohang15) + '\t' + str(FrequencyOf_allele1_Geoje15) + '\t' + str(FrequencyOf_allele1_Namhae15) + '\t' + str(FrequencyOf_allele1_YellowSea16) + '\t' + str(FrequencyOf_allele1_Jukbyeon07) + '\t' + str(FrequencyOf_allele1_JinhaeBay07) + '\t' + str(FrequencyOf_allele1_JinhaeBay08) + '\t' + str(FrequencyOf_allele1_Boryeong07) + '\t' + str(FrequencyOf_allele1_Geoje14) + '\t' + str(FrequencyOf_allele1_Kodiak03) + '\t' + str(FrequencyOf_allele1_Adak06) + '\t' + str(FrequencyOf_allele1_WashCoast05) + '\t' + str(FrequencyOf_allele1_HecStrait04) + '\t' + str(FrequencyOf_allele1_SalishSea12_13) + '\t' + str(FrequencyOf_allele1_JuandeFuca12) + '\t' + str(FrequencyOf_allele1_PWSound12) + '\t' + str(FrequencyOf_allele1_UnimakPass03))
			bad_locus_write = str(bad_locus_freqs).replace('[','').replace(',','\t').replace(']', '').replace("'", '').replace(' ','').replace('\\n','').replace('\\t','\t')

			print bad_locus_write
			blacklisted_MAF.write(bad_locus_write + '\n')
			blacklisted_genotypes.write(mystring)

#close open files
genotypes_file.close()
blacklisted_genotypes.close()
blacklisted_MAF.close()
filtered_genotypes.close()
output_freqs.close()