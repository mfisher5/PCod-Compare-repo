#!/bin/bash

###### THIS SHELL SCRIPT WILL PRODUCE THE # OF SEQUENCES IN ALL SAM FILES IN A SPECIFIED FOLDER ######

##MF 8/23/2017

echo "You must have `samtools` installed to run this script. If you do not have samtools installed, exit now."
echo "--"
echo "What is the relative path to your .sam files?"
read PATH

cd $PATH

echo "finding all 'sam' files"
sam_file_array="$(find . -name '*.sam')"

for file in $sam_file_array
do
	echo $file >> sam_readcounts_b7.txt
done


for file in $sam_file_array
do
	samtools view -F $file >> sam_readcounts_b7.txt
done