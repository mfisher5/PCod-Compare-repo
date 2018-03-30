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
python 1/prep_for_extraFilters.py /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Compare-repo/scripts/PostStacksFiltering/batch_/mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Compare-repo/stacks_b8_wgenome_r05.catalog.snps.tsv /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Compare-repo/scripts/PostStacksFiltering/batch_/mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Compare-repo/stacks_b8_wgenome_r05.catalog.snps2.tsv /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Compare-repo/scripts/PostStacksFiltering/batch_/mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Compare-repo/stacks_b8_wgenome_r05.haplotypes.tsv /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Compare-repo/scripts/PostStacksFiltering/batch_/mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Compare-repo/stacks_b8_wgenome_r05.haplotypes2.tsv


#Preparing file for correcting genotypes
python /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Compare-repo/scripts/PostStacksFiltering/preparing_file_for_correcting_genotypes.py /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Compare-repo/scripts/PostStacksFiltering/batch_/mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Compare-repo/stacks_b8_wgenome_r05.haplotypes2.tsv /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Compare-repo/scripts/PostStacksFiltering/batch_/mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Compare-repo/stacks_b8_wgenome_r05.biallelic_catalog.tsv /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Compare-repo/scripts/PostStacksFiltering/batch_/mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Compare-repo/stacks_b8_wgenome_r05.catalog.snps2.tsv no

echo 'Did this produce the correct output?'
read ANSWER
if [ $ANSWER == 'no' ]; then
	exit 1
fi


#unzip all .tags and .matches files that are still gzipped

cd /mnt/hgfs/Pacific\ cod/DataAnalysis/PCod-Compare-repo/scripts/PostStacksFiltering

echo 'you are now in folder:'
pwd

