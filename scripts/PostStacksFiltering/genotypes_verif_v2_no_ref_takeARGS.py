############## EDITED FROM Marine Brieuc   mbrieuc <a> uw.edu -- TO TAKE COMMAND LINE ARGUMENTS

#MF 2/13/2017

# argument 1- absolute path to folder containing stacks files (DON'T INCLUDE FINAL "/")
# argument 2- population map file, used in stacks' populations
# argument 3- batch #
# argument 4- minimum stack depth
##########################################################################################
#python genotypes_verif_v2_no_ref.py

import sys
stacksDIRECT = sys.argv[1]
popmap = open(sys.argv[2], "r")
batch = sys.argv[3]
depth = str(int(sys.argv[4]) - 1)

## Open the output genotype file for writing
file_name = stacksDIRECT + "/batch_" + batch + ".CorrectedGenotypes_biallelic_depth" + depth + ".txt"
genotype_file=open(file_name ,'w')

## Catalog file with name of locus in first columns, haplotype1 and haplotype2 in second column, SNP positions in the following columns (prepared in prior script)
cat_file_name = stacksDIRECT + "/batch_" + batch + ".biallelic_catalog.tsv"
Catalog_file=open(cat_file_name,'r')

## Initiate dictionary: locus :: haplotype1,haplotype2
loci={}
for line in Catalog_file:						
	line_list=line.split('\t')
	locus=line_list[0]							## The locus name is the zeroth element of the list
	loci[locus] = line_list[1:]					## The SNPs are the first to xth elements of the list

## Create a list of loci that will be used later and sort the list 	
list_loci=loci.keys()
list_loci.sort()
#print list_loci 

												## Add a tab at the beginning of the output file and write all the loci (in order) on the first row, separated by a tab
genotype_file.write('\t')
for locus in list_loci:
	genotype_file.write("%s\t" %locus)
genotype_file.write('\n')


### THIS CODE SAVES EACH SHARED LOCUS (BETWEEN INDIVIDUALS AND THE CATALOG) AS A DICTIONARY KEY



# create the file with the list of individuals: 
indiv_file_name = stacksDIRECT + "/CorrectedGenos_Individs.txt"
list_of_indiv = open(indiv_file_name, "w")
for line in popmap: 
	linelist = line.strip().split()
	list_of_indiv.write(linelist[0] + "\n")
popmap.close()
list_of_indiv.close()

## Open the file with the list of the individuals
list_of_indiv=open(indiv_file_name,'r')
for indiv in list_of_indiv:						## For each string in the file
	indiv_name=indiv.split('\n')[0] 			# Split the string by the newline. Now you have a list called indiv_name. If you get an error message about this line, check your end of line characters in the individual_list file to make sure that they match the .split delimiter
	print indiv_name
	### Open the individual.tags.tsv file
	matches_file=open(stacksDIRECT + '/' + indiv_name +'.matches.tsv','r')	# Open the matches.tags file for each individual in indiv_name list
	ind_loc_vs_cat={}	

# Start a dictionary that will hold some info????
	cat_loc=[]
	next(matches_file)
	for line in matches_file:	          # For each line(string) in the matches.tsv file for that individual
		line_list=line.split('\t')				# Split that line by tabs. It is now a list. 
		if not line_list[2] in cat_loc:				# If the locus name in the second position of the line_list is not also in the catalog loci
			ind_num=line_list[4]				# The ID of this locus within this individual is stored in 4th element of list
			ind_loc_vs_cat[ind_num]=line_list[2]		# The ID of the catalog locus is in the second element of line_list
			cat_loc.append(line_list[2])			# Append this name to the cat_loc list	
	list_ind_loci=ind_loc_vs_cat.keys()	

# The list_in_loci become the keys to the dictionary you initialized
	#ind_loc_vs_cat = [int(x) for x in ind_loc_vs_cat]
	#ind_loc_vs_cat.sort()
#	print ind_loc_vs_cat
	matches_file.close()
	#print ind_loc_vs_cat['1']
	Indiv_file=open(stacksDIRECT + '/' + indiv_name +'.tags.tsv','r')	

	line_number=0	
	count_haplotype1=0
	count_haplotype2=0
	count_other=0
	keep_locus=0
	locus_name=""
	
	### create a few lists
	list_to_export=list()
	loci_found=list()
		
		
	### THIS CODE CORRECTS THE GENOTYPES. TO CALL A HETEROZYGOTE, EACH ALLELE MUST BE PRESENT TWICE AND THE TOTAL READ DEPTH OF THAT LOCUS MUST BE TEN. 
	### for each line in he individual.tags.tsv file
	for line in Indiv_file:
		line_number+=1		
		if "consensus" in line:	
			### if line has consensus in it, export information about the previous locus (because you know that you're done with the previous locus now that you've encountered a new locus_	
			if keep_locus ==1:
				### call genotypes, based on depth ### you can change the depth
				if (count_haplotype1 > 1 and count_haplotype2 > 1 and count_haplotype1+count_haplotype2>depth):
					genotype="ab"
				elif (count_haplotype2 > 1 and count_haplotype1+count_haplotype2>depth):
					genotype="b"
				elif (count_haplotype1 > 1 and count_haplotype1+count_haplotype2>depth):
					genotype="a"
				else:
					genotype="-"				
				#print genotype
				#print count_haplotype1
				#print count_haplotype2
				list_to_export.append([locus_name,indiv_name,genotype])
				#print list_to_export
			### Reinitialize the parameters
			count_haplotype1=0
			count_haplotype2=0
			count_other=0			
			line_list=line.split('\t')
			#print line_list
			loc_name_ind=line_list[2]
			#print line_list[2]
			if(loc_name_ind in list_ind_loci):
				# If locus in the catalog, extract haplotype and SNP positions
				if (ind_loc_vs_cat[loc_name_ind]) in list_loci:
					locus_name=ind_loc_vs_cat[loc_name_ind]
					loci_found.append(locus_name)
					keep_locus = 1
					locus_line = line_number
					haplotype1=loci[locus_name][0]
					#print haplotype1
					haplotype2=loci[locus_name][1]
					#print haplotype2
					position= list()
					RANGE = range(2,len(haplotype1)+2,1)
					for index in RANGE:
						position.append(loci[locus_name][index]) 

				else:
					#print (ind_loc_vs_cat[loc_name_ind] + ' not in list! oups!'	)					
					keep_locus=0
			else:
				keep_locus=0
		### If line doesn't have "consensus in it, look for haplotypes in the reads if it is a locus that is in the catalog
		elif (keep_locus==1 and line_number>(locus_line+1)):
			haplotype=""
			line_list=line.split('\t')			
			sequence=line_list[9]
			#print position[0]
			for where in position:
				if len(where)>0: # added this if statement because there were some empty sections DD 4/24/13
#				print line_list
#					print '"' + where + '"'
#				print 'seq=', sequence
#				print 'seq[where]=', sequence[int(where)]
#				print line
					haplotype=haplotype+sequence[int(where.split()[0])] # added the split stuff DD 4/24/13
			#print haplotype
			if haplotype==haplotype1:
				count_haplotype1 += 1
			elif haplotype==haplotype2:
				count_haplotype2 += 1
			else:
				count_other += 1
	
	### When at the end of the file, call the genotype if the last locus of the file was in the catalog
	if keep_locus ==1:
		if (count_haplotype1 > 1 and count_haplotype2 > 1 and count_haplotype1+count_haplotype2>depth):
			genotype="ab"
		elif (count_haplotype2 > 1 and count_haplotype1+count_haplotype2>depth):
			genotype="b"
		elif (count_haplotype1 > 1 and count_haplotype1+count_haplotype2>depth):
			genotype="a"
		else:
			genotype="-"					
		list_to_export.append([locus_name,indiv_name,genotype])
	
	### Identify all loci that were missing in the individual and write them as missing value
#	print list_to_export	
	loci_missing = list(set(list_loci)-set(loci_found))
	#print loci_missing
	for locus in loci_missing:
		list_to_export.append([locus,indiv_name,"-"])	
	### At this point, sort by locus, all the loci in the dictionary should also be in this list
	list_to_export.sort()
	#print list_to_export
	### If the list of loci in the dictionary = list of loci to export (it should, but it's just to double check), write the genotypes in the output file (each line = indiv_name followed by the genotypes
	loci_in_order=list()
	genotype_in_order=""
	for locus in list_to_export:
		loci_in_order.append(locus[0])
		genotype_in_order=genotype_in_order + locus[2] + '\t'
	if loci_in_order==list_loci:
		genotype_file.write(indiv_name + '\t' + genotype_in_order + '\n')	
			
	Indiv_file.close()
genotype_file.close()
list_of_indiv.close()
Catalog_file.close()
	





						
