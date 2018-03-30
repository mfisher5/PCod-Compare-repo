### Argument 1 : haplotype_file
### Argument 2 : output_file
### Argument 3 : catalog.snp.tsv file
### Argument 4 : minimum coverage required


import commands
import sys
import glob
import os


infile=open(sys.argv[1],"r")
outfile=open(sys.argv[2],"w")
position_infile=open(sys.argv[3],"r")

no_consensus=[]

for line in infile:
 	split_line=line.split('\t')
	#print split_line
	#print split_line[1]
	if (split_line[1] != '0' and split_line.count('consensus')==0):
		no_consensus.append(line)

'''for item in no_consensus:
	if not "consensus" in line:
		no_consensus.append(line)
'''

locus_positions=[]
loci_positions=[]
previous_name='0'
for line in position_infile:
	split_line=line.strip().split('\t')
	locus_name=split_line[2]
	position=split_line[3]
	if(locus_name==previous_name):
		locus_positions.append(position)
	else:
		loci_positions.append(locus_positions)
		locus_positions=[]
		locus_positions.append(locus_name)
		locus_positions.append(position)
	previous_name=locus_name
loci_positions.append(locus_positions)
#print loci_positions
###Initiate dictionary: locus :: position1, position2, position3, etc.
loci={}
number_loci_cat=len(loci_positions)
for i in (range(1,int(len(loci_positions))-1)):
	locus=loci_positions[i][0]
	loci[locus] = loci_positions[i][1:]
	#print loci[locus]



#print loci_positions[0]
#print loci_positions[2]
#print loci_positions[12169]
#print len(loci_positions)
#print (no_consensus[1]) 
no_consensus=[w.replace("/"," ") for w in no_consensus]
no_consensus=[w.replace("\t"," ") for w in no_consensus]
no_consensus=[w.replace("\n","") for w in no_consensus]
no_consensus=[w.replace(" -","") for w in no_consensus]
no_consensus=[w.replace("  ","") for w in no_consensus]
#print (no_consensus[1])

for i in range(1,20):
	print no_consensus[i]

all_haplotypes_no_consensus=[]
unique_haplotypes=[]
all_hap_to_keep=[]
loci_name_keep=[]
for item in no_consensus:
	list_items=item.split()
	#for haplotype in list_items:
		#print(haplotype)
	#print list_items
	name_locus=list_items[0]
	coverage=list_items[1]
	if coverage>=sys.argv[4]:	
		all_haplotypes_no_consensus.append(list_items)
		unique_hap = list(set(list_items[2:]))
		unique_hap.sort()
		number_alleles=len(unique_hap)
		#print number_alleles
		all_genotypes=list_items[2:]
		hap_to_keep=[]
		haplotype_to_keep=[]	
		for haplotype in (unique_hap):
			if(all_genotypes.count(haplotype)>1):
				hap_to_keep.append(haplotype)
				#print hap_to_keep 
		if(len(hap_to_keep)==2):
			haplotype_to_keep.append(name_locus)
			haplotype_to_keep.append(hap_to_keep)
			all_hap_to_keep.append(haplotype_to_keep)
			loci_name_keep.append(haplotype_to_keep[0])
			number_polym=len(hap_to_keep[0])
			hap_1=''
			hap_2=''
			good_positions=[]
			for i in range(0,number_polym):
				if(hap_to_keep[0][i]!=hap_to_keep[1][i]):
					hap_1=hap_1+hap_to_keep[0][i]
					hap_2=hap_2+hap_to_keep[1][i]
					good_positions.append(loci[name_locus][i])
			#print haplotype_to_keep[1] + " " + haplotype_to_keep[2] + "\n" + haplotype_to_keep[0] + " " + hap_1 + " " + hap_2 + str(good_positions) + "\n"
			outfile.write(name_locus + "\t" + hap_1 + "\t" + hap_2 )
			for i in range(0,len(good_positions)):
				outfile.write('\t'+ str(good_positions[i]))
			outfile.write('\n')

#print all_haplotypes_no_consensus[1:20]
infile.close()
outfile.close()
position_infile.close()
