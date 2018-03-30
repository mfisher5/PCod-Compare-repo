#!/bin/bash
### This shell script will run all of the lab's extra filtering scripts after populations ###
## M.Fisher 2/12/2017


echo 'Before running running this script, please be sure that you have the following:'
echo ''
echo '1. A single folder containing (1) stacks populations output, and (2) batch catalog files'
echo ''
echo '2. These additional python scripts: (1) prep_for_extraFilters.py, (2) preparing_file_for_correcting_genotypes.py, (3) gzip_MBgenotypesverif_BASHshell.sh, (4) MB_genotypes_verif_v2_no_ref.py, (5) genepop_conversion_corrected.py, (6) transpose.py, (7) Eleni_filter_by_MinorAlleleFrequency (8) FilterLoci_by_MissingValues.py.'
echo ''
echo '3. The INPUT text file, UPDATED FOR YOUR DATA'
echo ''
echo 'Are you ready for the script to run?'
read ANSWER
if [ $ANSWER == 'no' ]; then
	exit 1
fi


HOME=$(pwd)


#Fix .haplotypes.tsv and .catalog.snps.tsv files for filtering
python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/prep_for_extraFilters.py /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.catalog.snps.tsv /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.catalog.snps2.tsv /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.haplotypes.tsv /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.haplotypes2.tsv


#Preparing file for correcting genotypes
python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/preparing_file_for_correcting_genotypes.py /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.haplotypes2.tsv /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.biallelic_catalog.tsv /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.catalog.snps2.tsv 1

echo 'Did this produce the correct output?'
read ANSWER
if [ $ANSWER == 'no' ]; then
	exit 1
fi


#unzip all .tags and .matches files that are still gzipped

cd /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5

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



#call Marine's script, edited to take arguments at the command line

echo 'verifying genotypes...'
python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/genotypes_verif_v2_no_ref_takeARGS.py /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5 /mnt/hgfs/Pacific\ cod/DataAnalysis/scripts/PopMap_L1L2stacks.txt 2 1

echo 'Did this produce the correct output?'
read ANSWER
if [ $ANSWER == 'no' ]; then
	exit 1
fi


#Run genepop_conversion_corrected.py

python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/Genepop_conversion_corrected.py /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.CorrectedGenotypes_biallelic.txt /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.CorrectedGenotypes_biallelic.genepop

echo 'Did this produce the correct output?'
read ANSWER
if [ $ANSWER == 'no' ]; then
	exit 1
fi


#Transpose the file for MAF calculation

echo 'checking genepop file for sample column header'

python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/add_sample_to_genepop.py /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5 2

echo 'running transpose.py'
python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/transpose.py /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.CorrectedGenos_biallelic_edit.genepop /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.CorrectedGenos_biallelic_TRANSPOSED.genepop


#Convert to .csv for MAF filtering

echo 'creating csv file...'
python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/convert_genepop_to_csv.py /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5 2

echo 'Did this produce the correct output?'
read ANSWER
if [ $ANSWER == 'no' ]; then
	exit 1
fi


#Run Eleni's MAF filtering script

echo 'running Eleni_filter_by_MinorAlleleFrequency_takeARGS.py'
python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/genMAFfiltering.py /mnt/hgfs/Pacific\ cod/DataAnalysis/scripts/PopMap_L1L2stacks.txt

python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/Eleni_filter_by_MinorAlleleFrequency_takeARGS.py /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.CorrectedGenos_biallelic_TRANSPOSED.csv /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.batch_1.MAFfiltering_outputFreqs /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.batch_1.filteredMAF /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.batch_1.MAFfiltering_BADgenotypes /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.batch_1.MAFfiltering_blacklistedMAF

echo 'Did this produce the correct output?'
read ANSWER
if [ $ANSWER == 'no' ]; then
	exit 1
fi


#Filter remaining loci for Missing Values

echo 'running FilterLoci_by_MissingValues_takeARGS.py'
python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/genMissingLoci.py /mnt/hgfs/Pacific\ cod/DataAnalysis/scripts/PopMap_L1L2stacks.txt

python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/FilterLoci_by_MissingValues.py /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.filteredMAF /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.filteredMAF_filteredLoci /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.LociFiltering_blacklistedLoci

echo 'Congratulations! You have successfully filtered your final genepop file. To view your final file, go to: batch_X.filteredMAF_filteredLoci'
