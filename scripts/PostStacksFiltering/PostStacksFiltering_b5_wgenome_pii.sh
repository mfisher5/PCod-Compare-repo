#!/bin/bash

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
python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/genMAFfiltering.py /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PopMap_L1-4_mdFilter.txt

python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/Eleni_filter_by_MinorAlleleFrequency_takeARGS.py /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.CorrectedGenotypes_biallelic_TRANSPOSED.csv /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.MAFfiltering_outputFreqs /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.filteredMAF /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.MAFfiltering_BADgenotypes /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.MAFfiltering_blacklistedMAF

echo 'Did this produce the correct output? (yes/no)'
read ANSWER
if [ $ANSWER == 'no' ]; then
	exit 1
fi


#Filter remaining loci for Missing Values

echo 'running FilterLoci_by_MissingValues_takeARGS.py'
python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/genMissingLoci.py /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PopMap_L1-4_mdFilter.txt

python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/FilterLoci_by_MissingValues.py /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.filteredMAF /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.filteredMAF_filteredLoci /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/stacks_b5_wgenome/batch_5.LociFiltering_blacklistedLoci

echo 'Congratulations! You have successfully filtered your final genepop file. To view your final file, go to: batch_X.filteredMAF_filteredLoci'
