#!/bin/bash

###### THIS SHELL SCRIPT WILL PRODUCE THE # OF MAPPED SEQUENCES IN ALL SAM FILES IN A SPECIFIED FOLDER ######

##MF 8/23/2017

echo "You must have `samtools` installed to run this script. If you do not have samtools installed, exit now."
echo "--"
echo "What is the relative path to your .sam files?"
read DIRECTORY

cd $DIRECTORY

echo "finding all 'sam' files"
sam_file_array="$(find . -name '*.sam.gz')"


echo "counting all reads / mapped reads in each samfile..."

for file in $sam_file_array
do
	gzip -d $file
	unzippedfiles="$(find . -name '*.sam')"
	for newfile in $unzippedfiles
	do
		echo $newfile >> ../sam_mappedreadcounts_b7_KOR.txt
		samtools view $newfile | wc -l >> ../sam_mappedreadcounts_b7_KOR.txt
		samtools view -S -F4 $newfile | wc -l >> ../sam_mappedreadcounts_b7_KOR.txt
		gzip $newfile
	done
	echo $file "done"
done
