## This script was written by Dan Drinan and transposes the rows and columns of a haplotype file
#python transpose.py batch_2_corrected_genotypes_2_alleles_genepop.txt batch_2_corrected_genotypes_2_alleles_genepop_transposed.txt

import sys

input_file = open(sys.argv[1], 'r')
input_file.readline() # header

matrix_of_data = []

for line in input_file:
	tmp_line = ''
	tmp_line += line
	matrix_of_data.append(tmp_line.split(","))

input_file.close()
print len(matrix_of_data)


transposed = zip(*matrix_of_data)
print len(transposed)

output_file = open(sys.argv[2], 'w')


for line in transposed:
    tmp_output = str(line).replace('[', '').replace('(', '').replace("'", '').replace(')', '').replace(',', '').replace(']', '') + '\n'
    output_file.write(tmp_output)

output_file.close()


