

################################################################################
#
# calc_het.py - this is a script used to calculate heterozygosity, based off of Dan Drinan's hetVsReadDepth.py
#
# 2017-March-28
# Daniel Drinan (ddrinan@uw.edu)
# MF Edited 5/9/2017

################################################################################

import argparse, subprocess

parser = argparse.ArgumentParser()
# it is either run with '-l' if you want to compare heterozygosity and read depth
# of a bunch of individuals
#
# or, if you are only interested in a single individual, you run with '-i' and
# '-f' 
parser.add_argument("-l", "--list", help="Population map, or any white space delimited list of individuals with sample name in first column")
parser.add_argument("-i", "--ind", help="name of individual to investigate (mutually \
                                         exclusive to '-l' and requires '-f') - \
                                         UNTESTED")
parser.add_argument("-f", "--file", help="location of file with genotypes (assumes \
                                          genepop format)")
parser.add_argument("-o", "--output", help="name of output file")
parser.add_argument("-d", "--denominator", help="use 2 if counting a FASTA file \
                                                 or 4 if counting a FASTQ file")
args = parser.parse_args()

output_file = open(args.output, 'w')





########################

# function to count the proportion of heterozygous loci in a sample
def countHet(sample_name):
        individuals_genotypes = subprocess.Popen(["grep " + sample_name + " " + \
                                args.file], stdout=subprocess.PIPE, shell=True)
        (genotypes_out, genotypes_err) = individuals_genotypes.communicate()
        genotypes_out = genotypes_out.split(',')[1] # removing everything except genotypes
        genotypes_out = genotypes_out.split()

        tmp_het = 0.0 # number of heterozygotes
        tmp_total = 0.0000000000000001 # number of genotyped loci

        for item in genotypes_out:
            if item.count('0') < len(item): # if true, a genotype exists
                tmp_total += 1

                # split the genotype in half and compare each allele
                if item[0:len(item)/2] != item[len(item)/2:]: # if true, it is a het
                    tmp_het += 1

        return tmp_het/tmp_total



#############################
##
## main
##

#############################
output_file.write('sample prop_het\n')

if args.list:

    list_file = open(args.list, 'r')

    for line in list_file:
        sample_name = line.strip().split()[0]

        # extract the list of genotypes for the individual from the genepop file

        tmp_proportion_het = countHet(sample_name)

        tmp_output = sample_name + ' ' + str(tmp_proportion_het) + '\n'

        output_file.write(tmp_output)


    list_file.close()

output_file.close()