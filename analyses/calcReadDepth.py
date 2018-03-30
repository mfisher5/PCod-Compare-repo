

################################################################################
#
# hetVsReadDepth.py - this is a script used to compare heterozygosity vs.
# read depth in an individual or list of individuals
#
# 2017-March-28
# Daniel Drinan (ddrinan@uw.edu)
#
################################################################################

import argparse, subprocess

parser = argparse.ArgumentParser()
# it is either run with '-l' if you want to compare heterozygosity and read depth
# of a bunch of individuals
#
# or, if you are only interested in a single individual, you run with '-i' and
# '-f' 
parser.add_argument("-l", "--list", help="white space delimited list of individuals \
                                          (column1) and corresponding sequence file \
                                          (column2; FASTA or FASTQ) to investigate (must include relative path to file")
parser.add_argument("-i", "--ind", help="name of individual to investigate (mutually \
                                         exclusive to '-l' and requires '-f') - \
                                         UNTESTED")
parser.add_argument("-o", "--output", help="name of output file")
parser.add_argument("-d", "--denominator", help="use 2 if counting a FASTA file \
                                                 or 4 if counting a FASTQ file")
args = parser.parse_args()

output_file = open(args.output, 'w')





#############################
## 
## functions
##
#############################

# function to count the number of lines in a fastq file
# the function returns the number of sequences counted
def countFastq(filename):
    denominator = 4
    if args.denominator == '2':
        denominator = 2

    num_lines_fastq = subprocess.Popen(["cat " + filename + " | wc -l"], \
                                       stdout=subprocess.PIPE, shell=True)
    (fastq_out, fastq_err) = num_lines_fastq.communicate()
    fastq_out = float(fastq_out.split()[0])/denominator # divide by 2 or 4 depending 
                                                        # on file type

    return fastq_out



#############################
##
## main
##
#############################
output_file.write('sample num_seq\n')

if args.list:

    list_file = open(args.list, 'r')

    for line in list_file:
        sample_name, fastq_file = line.split()
        print 'Processing', sample_name

        # count the number of lines in the FASTQ file
        tmp_num_lines = countFastq(fastq_file)

        tmp_output = sample_name + ' ' + str(tmp_num_lines) + '\n'

        output_file.write(tmp_output)


    list_file.close()

output_file.close()


