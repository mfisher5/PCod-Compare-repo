######## COUNT NUMBER OF LOCI PER SAMPLE AFTER SSTACKS ######
#
# MF 1/15/2017
#
# This script will sort through your `.matches` files produced from sstacks and return the number of loci in each individual. 
# The script produces output files that can then be fed into R and plotted as a histogram by population (using Eleni's Scatterplot_loci_pstacks.R script)
#
# Note that there is an optional argument included. This argument can be used if you wish to see how the number of loci per individual changes depending on the minimum stack depth. The default for this optional argument is 10. 
#
##############################


import argparse

parser = argparse.ArgumentParser(description='Counting loci in individuals after sstacks')

parser.add_argument('-i', '--popmap', help="population map input file. will be used to split individuals by population")
parser.add_argument('-p', '--path', help="path (relative or absolute) to sstacks' .matches files. Do not include final '/' on path.")
parser.add_argument('-o', '--output', help="file that will be used to store loci counts and other output. DO NOT INCLUDE .TXT ON THE END of name.")
parser.add_argument('-d', '--depth', help="minimum stack depth required to count a locus in an individual", default = '10')

args = parser.parse_args()


## set up from optional argument `--depth`
min_depth = int(args.depth)
depth_for_boolean = min_depth - 1


## using the population map, create a list of samples and a dictionary where you can look up each sample's population

popmap = open(args.popmap, "r")
samples = []
pop_dict = {}
for line in popmap:
    samples.append(line.strip().split()[0])    # add to list of samples
    pop_dict[line.strip().split()[0]] = line.strip().split()[1]    # create new entry in dictionary; key = sample, value = population
popmap.close()


## count the number of loci in each individual's matches file. will count total loci and loci with depth > `-d`

n_loci_dict = {}
n_loci_min_depth_dict = {}
for sample in samples:
    matchesfile = open(args.path + "/" + sample + ".matches.tsv", "r")
    matchesfile.readline() #header
    sample_loci_dict = {}
    for line in matchesfile:
        linelist = line.strip().split()
        cat_locus = linelist[2]
        depth = int(linelist[6])
        if cat_locus not in sample_loci_dict.keys():
            sample_loci_dict[cat_locus] = depth
        elif cat_locus in sample_loci_dict.keys():
            dict_depth = sample_loci_dict[cat_locus]
            dict_depth += depth
            sample_loci_dict[cat_locus] = dict_depth
    matchesfile.close()
    n_loci_dict[sample] = len(sample_loci_dict.keys())
    n_loci_min_depth = 0
    for locus in sample_loci_dict.keys():
        if sample_loci_dict[locus] > depth_for_boolean:
            n_loci_min_depth += 1
    n_loci_min_depth_dict[sample] = n_loci_min_depth   


## write outfiles

### for all loci
outfile = open(args.output + ".txt", "w")

outfile.write("sample\tn_loci\tpopulation\n")

for sample in n_loci_dict.keys():
    population = pop_dict[sample]
    outfile.write(sample + "\t" + str(n_loci_dict[sample]) + "\t" + population + "\n")
outfile.close()


### just for loci with minimum stack depth 
outfile = open(args.output + "_mindepth.txt", "w")

outfile.write("sample\tn_loci\tpopulation\n")

for sample in n_loci_min__dict.keys():
    population = pop_dict[sample]"
    outfile.write(sample + "\t" + str(n_loci_min_depth_dict[sample]) + "\t" + population + "\n")
outfile.close()
      

