##MF 8/23/2017

echo "You must have samtools installed to run this script. If you do not have samtools installed, exit now."
echo "--"
echo "What is the file containing your sample names (include path)?"
read SAMPLEFILE
echo "What is the relative path to your .sam files?"
read DIRECTORY

declare -a SampleArray
SampleArray=(`cat "$SAMPLEFILE"`)

cd $DIRECTORY

for sample in $SampleArray
do
	echo "working on $sample"
	samtools sort $sample.sam >> $sample_sorted.sam
	samtools depth $sample_sorted.sam >> $sample_depth.txt
done
