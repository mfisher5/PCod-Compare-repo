#!/bin/bash

###### THIS SHELL SCRIPT WILL PRODUCE THE # OF MAPPED SEQUENCES IN ALL SAM FILES IN A SPECIFIED FOLDER ######

##MF 8/23/2017

echo "You must have `samtools` installed to run this script. If you do not have samtools installed, exit now."
echo "--"
echo "What is the relative path to your .sam files?"
read DIRECTORY

cd $DIRECTORY

echo "finding all 'sam' files"
sam_file_array="$(find . -name '*.sam')"


echo "counting all reads / mapped reads in each samfile..."

for file in $sam_file_array
do
	echo $file >> ../sam_mappedreadcounts_b2_KORsubset.txt
	samtools view $file | wc -l >> ../sam_mappedreadcounts_b2_KORsubset.txt
	samtools view -S -F4 $file | wc -l >> ../sam_mappedreadcounts_b2_KORsubset.txt
	echo $file "done"
done
