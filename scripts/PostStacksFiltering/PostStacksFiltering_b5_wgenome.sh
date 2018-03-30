#!/bin/bash
### This shell script will run all of the lab's extra filtering scripts after populations ###
## M.Fisher 2/12/2017


echo 'Before running running this script, please be sure that you have the following:'
echo ''
echo '1. A single folder containing (1) stacks populations output, and (2) batch catalog files'
echo ''
echo '2. All additional python scripts listed in the input text file.'
echo ''
echo '3. The INPUT text file, UPDATED FOR YOUR DATA'
echo ''
echo 'Are you ready for the script to run? (yes/no)'
read ANSWER
if [ $ANSWER == 'no' ]; then
	exit 1
fi


HOME= /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering


#Fix .haplotypes.tsv and .catalog.snps.tsv files for filtering
python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/prep_for_extraFilters.py /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.catalog.snps.tsv /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.catalog.snps2.tsv /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.haplotypes.tsv /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.haplotypes2.tsv


#Preparing file for correcting genotypes
python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/preparing_file_for_correcting_genotypes.py /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.haplotypes2.tsv /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.biallelic_catalog.tsv /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.catalog.snps2.tsv 1

echo 'Did this produce the correct output? (yes/no)'
read ANSWER
if [ $ANSWER == 'no' ]; then
	exit 1
fi


#unzip all .tags and .matches files that are still gzipped

cd /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome

echo 'you are now in folder:'
pwd

echo "finding all gzipped tags files"

tags_file_array="$(find . -name '*.tags.tsv.gz')"

echo '--'
echo 'unzipping tags files...'

for file in $tags_file_array
do
	echo $file
	gzip -d $file
	echo 'file unzipped'
done

echo '--'

matches_file_array="$(find . -name '*.matches.tsv.gz')"

echo 'unzipping matches files...'

for file in $matches_file_array
do
	echo $file
	gzip -d $file
	echo 'file unzipped'
done


cd $HOME

echo 'Returning to home directory...'


#call Marine's script, edited to take arguments at the command line

echo 'verifying genotypes...'
python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/genotypes_verif_v2_no_ref_takeARGS.py /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PopMap_L1-4.txt 5 1

echo 'Did this produce the correct output? (yes/no)'
read ANSWER
if [ $ANSWER == 'no' ]; then
	exit 1
fi


#Run genepop_conversion_corrected.py

python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/Genepop_conversion_corrected.py /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.CorrectedGenotypes_biallelic.txt /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.CorrectedGenotypes_biallelic.genepop

echo 'Did this produce the correct output? (yes/no)'
read ANSWER
if [ $ANSWER == 'no' ]; then
	exit 1
fi

echo 'Before moving on to the MAF and loci missing data filtering, YOU MAY NEED TO FILTER FOR INDIVIDUALS MISSING DATA AT THIS TIME.'
echo 'Use the file batch_5.CorrectedGenotypes_biallelic.genepop to filter out individuals.'
echo 'Create a new population map file containing only the retained individuals.'
echo 'Then, when you are ready to continue, type in the filtered file name (without directory; for example, batch_1.CorrectedGenotypes_biallelic_filtered.genepop)'
read FILENAME
#Transpose the file for MAF calculation

echo 'Preparing for MAF filtering'

echo 'checking your genepop file for sample column header'

python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/add_sample_to_genepop_filteredindivids.py /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome 5 $FILENAME

echo 'running transpose.py'
python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/transpose.py /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.CorrectedGenotypes_biallelic_edit.genepop /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.CorrectedGenotypes_biallelic_TRANSPOSED.genepop


#Convert to .csv for MAF filtering

echo 'creating csv file...'
python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/convert_genepop_to_csv.py /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome 5

echo 'Did this produce the correct output? (yes/no)'
read ANSWER
if [ $ANSWER == 'no' ]; then
	exit 1
fi


#Run Eleni's MAF filtering script

echo 'running Eleni_filter_by_MinorAlleleFrequency_takeARGS.py'
python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/genMAFfiltering.py /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PopMap_L1-4_mdFilter

python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/Eleni_filter_by_MinorAlleleFrequency_takeARGS.py /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.CorrectedGenotypes_biallelic_TRANSPOSED.csv /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.MAFfiltering_outputFreqs /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.filteredMAF /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.MAFfiltering_BADgenotypes /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.MAFfiltering_blacklistedMAF

echo 'Did this produce the correct output? (yes/no)'
read ANSWER
if [ $ANSWER == 'no' ]; then
	exit 1
fi


#Filter remaining loci for Missing Values

echo 'running FilterLoci_by_MissingValues_takeARGS.py'
python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/genMissingLoci.py /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PopMap_L1-4_mdFilter

python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/FilterLoci_by_MissingValues.py /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.filteredMAF /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.filteredMAF_filteredLoci /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.LociFiltering_blacklistedLoci

echo 'Congratulations! You have successfully filtered your final genepop file. To view your final file, go to: batch_X.filteredMAF_filteredLoci'
