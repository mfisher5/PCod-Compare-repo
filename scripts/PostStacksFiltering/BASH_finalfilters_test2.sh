#!/bin/bash
### This shell script will run all of the lab's extra filtering scripts after populations ###
## M.Fisher 2/12/2017

#Run Eleni's MAF filtering script

echo 'running Eleni_filter_by_MinorAlleleFrequency_takeARGS.py'
python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/genMAFfiltering.py /mnt/hgfs/Pacific\ cod/DataAnalysis/scripts/PopMap_L1L2stacks.txt

python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Korea-repo/scripts/PostStacksFiltering/Eleni_filter_by_MinorAlleleFrequency_takeARGS.py /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.CorrectedGenotypes_biallelic_TRANSPOSED.csv /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.MAFfiltering_outputFreqs /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.filteredMAF /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.MAFfiltering_BADgenotypes /mnt/hgfs/Pacific\ cod/DataAnalysis/L1L2stacks_m5/batch_2.MAFfiltering_blacklistedMAF

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
