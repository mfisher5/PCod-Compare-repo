## Sort a .sam file, and then determine the coverage at each catalog locus
##MF 9/6/2017

echo "You must have samtools installed to run this script. If you do not have samtools installed, exit now."
echo "--"
echo "What is the relative path to your .sam files?"
read DIRECTORY

declare -a SampleArray=(KOD03_035 KOD03_051 KOD03_052 KOD03_053 KOD03_054 AD06_001 AD06_002 AD06_003 AD06_004 AD06_005 WC05_001 WC05_002 WC05_003 WC05_004 WC05_005 HS04_001 HS04_002 HS04_003 HS04_004 HS04_005 PS12_001 PS12_002 PS12_003 PS12_004 PS12_005 GS13_001 GS13_002 GS13_003 GS13_004 GS13_005 PWS12_001 PWS12_002 PWS12_003 PWS12_004 PWS12_005 UP03_001 UP03_009 UP03_017 UP03_025 UP03_033)


cd $DIRECTORY

for sample in $SampleArray
do
	echo "working on $sample"
	samtools sort $sample.sam >> $sample-sorted.sam
	samtools depth $sample-sorted.sam >> $sample-depth.txt
done
