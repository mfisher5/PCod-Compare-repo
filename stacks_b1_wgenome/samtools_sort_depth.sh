#!/bin/bash

## Sort a .sam file, and then determine the coverage at each catalog locus
#MF 9/6/2017

samtools sort KOD03_035.sam >> KOD03_035_sorted.sam
samtools depth KOD03_035_sorted.sam >> KOD03_035_depth.txt
rm KOD03_035_sorted.sam
echo 'KOD03_035done.'

samtools sort KOD03_051.sam >> KOD03_051_sorted.sam
samtools depth KOD03_051_sorted.sam >> KOD03_051_depth.txt
rm KOD03_051_sorted.sam
echo 'KOD03_051done.'

samtools sort KOD03_052.sam >> KOD03_052_sorted.sam
samtools depth KOD03_052_sorted.sam >> KOD03_052_depth.txt
rm KOD03_052_sorted.sam
echo 'KOD03_052done.'

samtools sort KOD03_053.sam >> KOD03_053_sorted.sam
samtools depth KOD03_053_sorted.sam >> KOD03_053_depth.txt
rm KOD03_053_sorted.sam
echo 'KOD03_053done.'

samtools sort KOD03_054.sam >> KOD03_054_sorted.sam
samtools depth KOD03_054_sorted.sam >> KOD03_054_depth.txt
rm KOD03_054_sorted.sam
echo 'KOD03_054done.'

samtools sort AD06_001.sam >> AD06_001_sorted.sam
samtools depth AD06_001_sorted.sam >> AD06_001_depth.txt
rm AD06_001_sorted.sam
echo 'AD06_001done.'

samtools sort AD06_002.sam >> AD06_002_sorted.sam
samtools depth AD06_002_sorted.sam >> AD06_002_depth.txt
rm AD06_002_sorted.sam
echo 'AD06_002done.'

samtools sort AD06_003.sam >> AD06_003_sorted.sam
samtools depth AD06_003_sorted.sam >> AD06_003_depth.txt
rm AD06_003_sorted.sam
echo 'AD06_003done.'

samtools sort AD06_004.sam >> AD06_004_sorted.sam
samtools depth AD06_004_sorted.sam >> AD06_004_depth.txt
rm AD06_004_sorted.sam
echo 'AD06_004done.'

samtools sort AD06_005.sam >> AD06_005_sorted.sam
samtools depth AD06_005_sorted.sam >> AD06_005_depth.txt
rm AD06_005_sorted.sam
echo 'AD06_005done.'

samtools sort WC05_001.sam >> WC05_001_sorted.sam
samtools depth WC05_001_sorted.sam >> WC05_001_depth.txt
rm WC05_001_sorted.sam
echo 'WC05_001done.'

samtools sort WC05_002.sam >> WC05_002_sorted.sam
samtools depth WC05_002_sorted.sam >> WC05_002_depth.txt
rm WC05_002_sorted.sam
echo 'WC05_002done.'

samtools sort WC05_003.sam >> WC05_003_sorted.sam
samtools depth WC05_003_sorted.sam >> WC05_003_depth.txt
rm WC05_003_sorted.sam
echo 'WC05_003done.'

samtools sort WC05_004.sam >> WC05_004_sorted.sam
samtools depth WC05_004_sorted.sam >> WC05_004_depth.txt
rm WC05_004_sorted.sam
echo 'WC05_004done.'

samtools sort WC05_005.sam >> WC05_005_sorted.sam
samtools depth WC05_005_sorted.sam >> WC05_005_depth.txt
rm WC05_005_sorted.sam
echo 'WC05_005done.'

samtools sort HS04_001.sam >> HS04_001_sorted.sam
samtools depth HS04_001_sorted.sam >> HS04_001_depth.txt
rm HS04_001_sorted.sam
echo 'HS04_001done.'

samtools sort HS04_002.sam >> HS04_002_sorted.sam
samtools depth HS04_002_sorted.sam >> HS04_002_depth.txt
rm HS04_002_sorted.sam
echo 'HS04_002done.'

samtools sort HS04_003.sam >> HS04_003_sorted.sam
samtools depth HS04_003_sorted.sam >> HS04_003_depth.txt
rm HS04_003_sorted.sam
echo 'HS04_003done.'

samtools sort HS04_004.sam >> HS04_004_sorted.sam
samtools depth HS04_004_sorted.sam >> HS04_004_depth.txt
rm HS04_004_sorted.sam
echo 'HS04_004done.'

samtools sort HS04_005.sam >> HS04_005_sorted.sam
samtools depth HS04_005_sorted.sam >> HS04_005_depth.txt
rm HS04_005_sorted.sam
echo 'HS04_005done.'

samtools sort PS12_001.sam >> PS12_001_sorted.sam
samtools depth PS12_001_sorted.sam >> PS12_001_depth.txt
rm PS12_001_sorted.sam
echo 'PS12_001done.'

samtools sort PS12_002.sam >> PS12_002_sorted.sam
samtools depth PS12_002_sorted.sam >> PS12_002_depth.txt
rm PS12_002_sorted.sam
echo 'PS12_002done.'

samtools sort PS12_003.sam >> PS12_003_sorted.sam
samtools depth PS12_003_sorted.sam >> PS12_003_depth.txt
rm PS12_003_sorted.sam
echo 'PS12_003done.'

samtools sort PS12_004.sam >> PS12_004_sorted.sam
samtools depth PS12_004_sorted.sam >> PS12_004_depth.txt
rm PS12_004_sorted.sam
echo 'PS12_004done.'

samtools sort PS12_005.sam >> PS12_005_sorted.sam
samtools depth PS12_005_sorted.sam >> PS12_005_depth.txt
rm PS12_005_sorted.sam
echo 'PS12_005done.'

samtools sort GS13_001.sam >> GS13_001_sorted.sam
samtools depth GS13_001_sorted.sam >> GS13_001_depth.txt
rm GS13_001_sorted.sam
echo 'GS13_001done.'

samtools sort GS13_002.sam >> GS13_002_sorted.sam
samtools depth GS13_002_sorted.sam >> GS13_002_depth.txt
rm GS13_002_sorted.sam
echo 'GS13_002done.'

samtools sort GS13_003.sam >> GS13_003_sorted.sam
samtools depth GS13_003_sorted.sam >> GS13_003_depth.txt
rm GS13_003_sorted.sam
echo 'GS13_003done.'

samtools sort GS13_004.sam >> GS13_004_sorted.sam
samtools depth GS13_004_sorted.sam >> GS13_004_depth.txt
rm GS13_004_sorted.sam
echo 'GS13_004done.'

samtools sort GS13_005.sam >> GS13_005_sorted.sam
samtools depth GS13_005_sorted.sam >> GS13_005_depth.txt
rm GS13_005_sorted.sam
echo 'GS13_005done.'

samtools sort PWS12_001.sam >> PWS12_001_sorted.sam
samtools depth PWS12_001_sorted.sam >> PWS12_001_depth.txt
rm PWS12_001_sorted.sam
echo 'PWS12_001done.'

samtools sort PWS12_002.sam >> PWS12_002_sorted.sam
samtools depth PWS12_002_sorted.sam >> PWS12_002_depth.txt
rm PWS12_002_sorted.sam
echo 'PWS12_002done.'

samtools sort PWS12_003.sam >> PWS12_003_sorted.sam
samtools depth PWS12_003_sorted.sam >> PWS12_003_depth.txt
rm PWS12_003_sorted.sam
echo 'PWS12_003done.'

samtools sort PWS12_004.sam >> PWS12_004_sorted.sam
samtools depth PWS12_004_sorted.sam >> PWS12_004_depth.txt
rm PWS12_004_sorted.sam
echo 'PWS12_004done.'

samtools sort PWS12_005.sam >> PWS12_005_sorted.sam
samtools depth PWS12_005_sorted.sam >> PWS12_005_depth.txt
rm PWS12_005_sorted.sam
echo 'PWS12_005done.'

samtools sort UP03_001.sam >> UP03_001_sorted.sam
samtools depth UP03_001_sorted.sam >> UP03_001_depth.txt
rm UP03_001_sorted.sam
echo 'UP03_001done.'

samtools sort UP03_009.sam >> UP03_009_sorted.sam
samtools depth UP03_009_sorted.sam >> UP03_009_depth.txt
rm UP03_009_sorted.sam
echo 'UP03_009done.'

samtools sort UP03_017.sam >> UP03_017_sorted.sam
samtools depth UP03_017_sorted.sam >> UP03_017_depth.txt
rm UP03_017_sorted.sam
echo 'UP03_017done.'

samtools sort UP03_025.sam >> UP03_025_sorted.sam
samtools depth UP03_025_sorted.sam >> UP03_025_depth.txt
rm UP03_025_sorted.sam
echo 'UP03_025done.'

samtools sort UP03_033.sam >> UP03_033_sorted.sam
samtools depth UP03_033_sorted.sam >> UP03_033_depth.txt
rm UP03_033_sorted.sam
echo 'UP03_033done.'

