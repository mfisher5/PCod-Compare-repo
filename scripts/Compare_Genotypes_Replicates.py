################## Code to compare replicates or reruns of the same sample ##################
#
# MF 4/23/2018
# built for comparison of subsetted / full data samples for simulation of low stack depth
#
##############################################################################################

import argparse 

parser = argparse.ArgumentParser(description="compare genotypes between an original sample and its replicate. You must have two matrix files, one for original samples and one for replicates.")

parser.add_argument("-i", "--input", help="first input file; genotypes in a 2 x 2 matrix; individuals by row, loci by columns")
parser.add_argument("-r", "--replicates", help="second input file with replicates; genotypes in a 2 x 2 matrix; individuals by row, loci by columns")
parser.add_argument("-by", "--by", help="do you want mismatches compared for each locus, or for each individual sample? [locus/sample]")
parser.add_argument("-o", "--outfile", help="prefix for output file")


args = parser.parse_args()


### open files ###
batch1 = open(args.input, "r")
batch2 = open(args.replicates, "r")
output = open(args.outfile, "w")

### call packages ###
import collections
import sys

########################### BY LOCUS ####################################

if args.by == "locus":
	print "You have chosen to calculate genotype error by locus."
	## initiate ordered dictionary so loci can be called with indices
	loci_dict_b1 = collections.OrderedDict()
	loci_dict_b2 = collections.OrderedDict()

	## read loci in as keys
	loci_list1 = batch1.readline().strip().split()
	loci_list2 = batch2.readline().strip().split()
	if len(loci_list1) != len(loci_list2):
		sys.exit("ERROR: You do not have the same number of loci in your files!")
	## fill in dictionary with genotypes (uses indices to assign genotypes in row to loci in dictionary)
	n_samples = 0
	for line in batch1:
    		genotypes = line.strip().split()[1:]
    		if n_samples == 0:
       			for i in range(0,len(genotypes)):
            			loci_dict_b1[loci_list1[i]] = [genotypes[i]]
        		n_samples += 1
    		elif n_samples > 0:
        		for i in range(0,len(genotypes)):
            			geno_list = loci_dict_b1[loci_list1[i]]
            			geno_list.append(genotypes[i])
            			loci_dict_b1[loci_list1[i]] = geno_list
        		n_samples += 1
	print "processed ", n_samples, " in first file."
	batch1.close()

	n_samples = 0
	for line in batch2:
    		genotypes = line.strip().split()[1:]
    		if n_samples == 0:
        		for i in range(0,len(genotypes)):
            			loci_dict_b2[loci_list2[i]] = [genotypes[i]]
        		n_samples += 1
    		elif n_samples > 0:
        		for i in range(0,len(genotypes)):
            			geno_list = loci_dict_b2[loci_list2[i]]
            			geno_list.append(genotypes[i])
           			loci_dict_b2[loci_list2[i]] = geno_list
       			n_samples += 1
	print "processed ", n_samples, " in second file."
	batch2.close()
	## compare genotypes at each locus
	geno_codes = {}
	for locus in loci_list1:
		b1_genotypes = loci_dict_b1[locus]
    		b2_genotypes = loci_dict_b2[locus]
    		coded = []
    		## add appropriate code to list based on comparison of two genotypes
    		#--- 0 = both missing
    		#--- 1 = batch 1 missing
    		#--- 2 = batch 2 missing
		#--- 3 = both genotyped
		#--- 4 = both genotyped, matched
    		#--- 5 = both genotyped, mismatch because of hom --> hom
    		#--- 6 = both genotyped, mismatched because of hom --> het
    		#--- 7 = both genotyped, mismatched because of het --> hom
    		for i in range(0,len(b1_genotypes)):
        		if b1_genotypes[i] == b2_genotypes[i]:
            		# if genotypes matched and neither are missing
            			if b1_genotypes[i] != "0000":
                			coded.append(4)
            			elif b1_genotypes[i] == "0000":
                			coded.append(0)
        		# if genotypes don't match
        		elif b1_genotypes[i] != b2_genotypes[i]:
            			# if they are mismatched because one is missing
            			if b1_genotypes[i] == "0000":
                			coded.append(1)
            			elif b2_genotypes[i] == "0000":
                			coded.append(2)
            		# if they are mismatched and neither are missing
            		else:
                		# b1 het (b2 must be hom)
                		if b1_genotypes[i][0:2] != b1_genotypes[i][2:]:
                    			coded.append(7)
               			# b1 hom
                		elif b1_genotypes[i][0:2] == b1_genotypes[i][2:]:
                    			# b2 hom
                    			if b2_genotypes[i][0:2] == b2_genotypes[i][2:]:
                        			coded.append(5)
                    			# b2 het
                    			elif b2_genotypes[i][0:2] != b2_genotypes[i][2:]:
                        			coded.append(6)
		geno_codes[locus] = coded
    		if len(coded) < 50:
        		print "processed only ", len(coded), " samples at locus."

	## record mismatches
	outfile.write("locus\tn.pairs\tn.both.miss\tn.b1.miss\tn.b2.miss\tn.both.genod\tn.matched\tn.mismatch.hom2hom\tn.mismatch.het2hom\tn.mismatch.hom2het\n")
	for locus in loci_list1:
    		coded = geno_codes[locus]
    		n_pairs = len(coded)
    		n_miss = len([i for i in coded if i == 0])
    		n_missb1 = len([i for i in coded if i == 1])
    		n_missb2 = len([i for i in coded if i == 2])
    		n_matched = len([i for i in coded if i == 4])
    		n_mismatch_homhom = len([i for i in coded if i == 5])
    		n_mismatch_homhet = len([i for i in coded if i == 6])
    		n_mismatch_hethom = len([i for i in coded if i == 7])
    		n_genod = n_matched + n_mismatch_homhom + n_mismatch_homhet + n_mismatch_hethom
    		outfile.write(locus + "\t" + str(n_pairs) + "\t" + str(n_miss) + "\t" + str(n_missb1) + "\t" + str(n_missb2) + "\t" + str(n_genod) + "\t")
    		outfile.write(str(n_matched) + "\t" + str(n_mismatch_homhom) + "\t" + str(n_mismatch_homhet) + "\t" + str(n_mismatch_hethom) + "\n")
	outfile.close()
	print "Wrote output file for mismatched genotypes by locus."


########################### BY INDIVIDUAL ####################################
elif args.by == "sample":
	print "You have chosen to calculate genotype error by sample."
	## initiate ordered dictionary so loci can be called with indices
	samples_dict_b1 = collections.OrderedDict()
	samples_dict_b2 = collections.OrderedDict()

	## read loci in as keys
	batch1.readline()
	batch2.readline()
	for line in batch1:
    		linelist = line.strip().split()
    		samples_dict_b1[linelist[0]] = linelist[1:]
	batch1.close()

	for line in batch2:
    		linelist = line.strip().split()
   		sampID = linelist[0].strip("_subset")
    		samples_dict_b2[sampID] = linelist[1:]
	batch2.close()
	if len([i for i in samples_dict_b1.keys() if i in samples_dict_b2.keys()]) < len(samples_dict_b1.keys()):
		sys.exit("You do not have the same number of samples in your files")
	## compare genotypes in each individual
	samples_geno_codes = {}

	for sample in samples_dict_b1.keys():
    		b1_genotypes = samples_dict_b1[sample]
    		b2_genotypes = samples_dict_b2[sample]
    		coded = []
    		## add appropriate code to list based on comparison of two genotypes
    		#--- 0 = both missing
    		#--- 1 = batch 1 missing
   		#--- 2 = batch 2 missing
    		#--- 3 = both genotyped
    		#--- 4 = both genotyped, matched
    		#--- 5 = both genotyped, mismatch because of hom --> hom
    		#--- 6 = both genotyped, mismatched because of hom --> het
    		#--- 7 = both genotyped, mismatched because of het --> hom
    		for i in range(0,len(b1_genotypes)):
        		if b1_genotypes[i] == b2_genotypes[i]:
            		# if genotypes matched and neither are missing
            			if b1_genotypes[i] != "0000":
                			coded.append(4)
            			elif b1_genotypes[i] == "0000":
                			coded.append(0)
        		# if genotypes don't match
        		elif b1_genotypes[i] != b2_genotypes[i]:
            			# if they are mismatched because one is missing
            			if b1_genotypes[i] == "0000":
                			coded.append(1)
            			elif b2_genotypes[i] == "0000":
                			coded.append(2)
            		# if they are mismatched and neither are missing
            		else:
                		# b1 het (b2 must be hom)
                		if b1_genotypes[i][0:2] != b1_genotypes[i][2:]:
                    			coded.append(7)
                		# b1 hom
                		elif b1_genotypes[i][0:2] == b1_genotypes[i][2:]:
                    			# b2 hom
                    			if b2_genotypes[i][0:2] == b2_genotypes[i][2:]:
                        			coded.append(5)
                    			# b2 het
                    			elif b2_genotypes[i][0:2] != b2_genotypes[i][2:]:
                        			coded.append(6)
    		samples_geno_codes[sample] = coded
    		if len(coded) < len(b1_genotypes):
        		print "Oh no! processed only ", len(coded), " out of ", len(b1_genotypes)," loci for sample ", sample

	## calculate mismatches	
	outfile.write("sample\tn.loci\tn.both.miss\tn.b1.miss\tn.b2.miss\tn.both.genod\tn.matched\tn.mismatch.hom2hom\tn.mismatch.het2hom\tn.mismatch.hom2het\n")

	for sample in samples_geno_codes.keys():
    		coded = samples_geno_codes[sample]
    		n_loci = len(coded)
    		n_miss = len([i for i in coded if i == 0])
    		n_missb1 = len([i for i in coded if i == 1])
    		n_missb2 = len([i for i in coded if i == 2])
    		n_matched = len([i for i in coded if i == 4])
    		n_mismatch_homhom = len([i for i in coded if i == 5])
    		n_mismatch_homhet = len([i for i in coded if i == 6])
    		n_mismatch_hethom = len([i for i in coded if i == 7])
    		n_genod = n_matched + n_mismatch_homhom + n_mismatch_homhet + n_mismatch_hethom
    		outfile.write(sample + "\t" + str(n_loci) + "\t" + str(n_miss) + "\t" + str(n_missb1) + "\t" + str(n_missb2) + "\t" + str(n_genod) + "\t")
    		outfile.write(str(n_matched) + "\t" + str(n_mismatch_homhom) + "\t" + str(n_mismatch_homhet) + "\t" + str(n_mismatch_hethom) + "\n")
	outfile.close()









