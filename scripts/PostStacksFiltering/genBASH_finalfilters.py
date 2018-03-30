### THIS SCRIPT GENERATES THE BASH SHELL THAT WILL RUN ALL OF THE EXTRA FILTERING STEPS AFTER POPULATIONS ###

## arguments: input file, pre-formatted
##MF 2/12/2017

###################################################################################################

import sys
input = open(sys.argv[1], "r")

#set variable names from input file
input.readline()
NewScript = input.readline().split("\t")[0]
stacksDIRECT = input.readline().split("\t")[0]
BATCH = input.readline().split("\t")[0]
SNPzipped = input.readline().split("\t")[0]
COVERAGE = input.readline().split("\t")[0]
prepDIRECT = input.readline().split("\t")[0]
prep2DIRECT = input.readline().split("\t")[0]
gunzipDIRECT = input.readline().split("\t")[0]
verifDIRECT = input.readline().split("\t")[0]
convertDIRECT = input.readline().split("\t")[0]
addDIRECT = input.readline().split("\t")[0]
transDIRECT = input.readline().split("\t")[0]
tocsvDIRECT = input.readline().split("\t")[0]
mafgenDIRECT = input.readline().split("\t")[0]
mvalsDIRECT = input.readline().split("\t")[0]
POPMAP = input.readline().split("\t")[0]



#create new script, initiate with intro text
shell = open(NewScript, "w")
shell.write("#!/bin/bash" + "\n")
shell.write("### This shell script will run all of the lab's extra filtering scripts after populations ###\n## M.Fisher 2/12/2017\n\n\n")

shell.write("echo 'Before running running this script, please be sure that you have the following:'\necho ''\necho '1. A single folder containing (1) stacks populations output, and (2) batch catalog files'\necho ''\necho '2. These additional python scripts: (1) prep_for_extraFilters.py, (2) preparing_file_for_correcting_genotypes.py, (3) gzip_MBgenotypesverif_BASHshell.sh, (4) MB_genotypes_verif_v2_no_ref.py, (5) genepop_conversion_corrected.py, (6) transpose.py, (7) Eleni_filter_by_MinorAlleleFrequency (8) FilterLoci_by_MissingValues.py.'\necho ''\necho '3. The INPUT text file, UPDATED FOR YOUR DATA'\necho ''\n")

shell.write("echo 'Are you ready for the script to run?'\nread ANSWER\nif [ $ANSWER == 'no' ]; then\n\texit 1\nfi\n\n")

shell.write("\nHOME=$(pwd)\n\n")


#prepare for the Marine's filtering scripts by altering haplotypes and catalog.snps files: 
shell.write("\n#Fix .haplotypes.tsv and .catalog.snps.tsv files for filtering\n")
if SNPzipped == "yes": 
	shell.write("gzip -d " + stacksDIRECT + "/batch_" + BATCH + ".catalog.snps.tsv.gz\n")
shell.write("python " + prepDIRECT + "/prep_for_extraFilters.py " + stacksDIRECT + "/batch_" + BATCH + ".catalog.snps.tsv " + stacksDIRECT + "/batch_" + BATCH + ".catalog.snps2.tsv " + stacksDIRECT + "/batch_" + BATCH + ".haplotypes.tsv " + stacksDIRECT + "/batch_" + BATCH + ".haplotypes2.tsv\n\n")

#Marine's prep script
shell.write("\n#Preparing file for correcting genotypes\n")
shell.write("python " + prep2DIRECT + "/preparing_file_for_correcting_genotypes.py " + stacksDIRECT + "/batch_" + BATCH + ".haplotypes2.tsv " + stacksDIRECT + "/batch_" + BATCH + ".biallelic_catalog.tsv " + stacksDIRECT + "/batch_" + BATCH + ".catalog.snps2.tsv " + COVERAGE + "\n\n") 

shell.write("echo 'Did this produce the correct output?'\nread ANSWER\nif [ $ANSWER == 'no' ]; then\n\texit 1\nfi\n\n")

#Marine's genotypes_verif_v2_no_ref_takeARGS (edited to take arguments at command line)
shell.write("\n#unzip all .tags and .matches files that are still gzipped\n\n")
shell.write("cd " + stacksDIRECT + "\n\n")
shell.write("echo 'you are now in folder:'\npwd\n\n")

gunzip = open(gunzipDIRECT + "/gunzip_string.txt", "r")
gunzip_str = gunzip.read()
gunzip.close()

shell.write(gunzip_str + "\n\n")

shell.write("cd $HOME\n\necho 'Returning to home directory...'\n\n")

shell.write("\n#call Marine's script, edited to take arguments at the command line\n\n")
shell.write("echo 'verifying genotypes...'\n")
shell.write("python " + verifDIRECT + "/genotypes_verif_v2_no_ref_takeARGS.py " + stacksDIRECT + " " + POPMAP + " " + BATCH + " " + COVERAGE + "\n\n")


shell.write("echo 'Did this produce the correct output?'\nread ANSWER\nif [ $ANSWER == 'no' ]; then\n\texit 1\nfi\n\n")


#Charlie's genepop_conversion_corrected.py script
shell.write("\n#Run genepop_conversion_corrected.py\n\n")
shell.write("python " + convertDIRECT + "/Genepop_conversion_corrected.py " + stacksDIRECT + "/batch_" + BATCH + ".CorrectedGenotypes_biallelic.txt " + stacksDIRECT + "/batch_" + BATCH + ".CorrectedGenotypes_biallelic.genepop\n\n")


shell.write("echo 'Did this produce the correct output?'\nread ANSWER\nif [ $ANSWER == 'no' ]; then\n\texit 1\nfi\n\n")


#Transpose the file for MAF calculation
shell.write("\n#Transpose the file for MAF calculation\n\n")
shell.write("echo 'checking genepop file for sample column header'\n\n")
### If necessary, add in the word "sample" as the first column header
shell.write("python " + addDIRECT + "/add_sample_to_genepop.py " + stacksDIRECT + " " + BATCH + "\n\n")
### Run transpose.py
shell.write("echo 'running transpose.py'\n")
shell.write("python " + transDIRECT + "/transpose.py " + stacksDIRECT + "/batch_" + BATCH + ".CorrectedGenotypes_biallelic_edit.genepop " + stacksDIRECT + "/batch_" + BATCH + ".CorrectedGenotypes_biallelic_TRANSPOSED.genepop\n\n")


#Convert to .csv to filter for MAF
shell.write("\n#Convert to .csv for MAF filtering\n\n")
shell.write("echo 'creating csv file...'\n")
shell.write("python " + convertDIRECT  + "/convert_genepop_to_csv.py " + stacksDIRECT + " " + BATCH + "\n\n")


shell.write("echo 'Did this produce the correct output?'\nread ANSWER\nif [ $ANSWER == 'no' ]; then\n\texit 1\nfi\n\n")


#Run MAF filtering
shell.write("\n#Run Eleni's MAF filtering script\n\n")
shell.write("echo 'running Eleni_filter_by_MinorAlleleFrequency_takeARGS.py'\n")
shell.write("python " + mafgenDIRECT  + "/genMAFfiltering.py " + POPMAP + "\n\n")
shell.write("python " + mafgenDIRECT + "/Eleni_filter_by_MinorAlleleFrequency_takeARGS.py " + stacksDIRECT + "/batch_" + BATCH + ".CorrectedGenotypes_biallelic_TRANSPOSED.csv " + stacksDIRECT + "/batch_" + BATCH + ".MAFfiltering_outputFreqs " + stacksDIRECT + "/batch_" + BATCH + ".filteredMAF " + stacksDIRECT + "/batch_" + BATCH + ".MAFfiltering_BADgenotypes " + stacksDIRECT + "/batch_" + BATCH + ".MAFfiltering_blacklistedMAF")


shell.write("\n\necho 'Did this produce the correct output?'\nread ANSWER\nif [ $ANSWER == 'no' ]; then\n\texit 1\nfi\n\n")


#Filter Loci with "x"% missing values
shell.write("\n#Filter remaining loci for Missing Values\n\n")
shell.write("echo 'running FilterLoci_by_MissingValues_takeARGS.py'\n")
shell.write("python " + mvalsDIRECT + "/genMissingLoci.py " + POPMAP + "\n\n")
shell.write("python " + mvalsDIRECT + "/FilterLoci_by_MissingValues.py " + stacksDIRECT + "/batch_" + BATCH + ".filteredMAF " + stacksDIRECT + "/batch_" + BATCH + ".filteredMAF_filteredLoci " + stacksDIRECT + "/batch_" + BATCH + ".LociFiltering_blacklistedLoci")

shell.write("\n\necho 'Congratulations! You have successfully filtered your final genepop file. To view your final file, go to: batch_X.filteredMAF_filteredLoci'")
