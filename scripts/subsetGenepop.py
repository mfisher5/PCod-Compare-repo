





################################################################################
#   
# subsetGenepop.py
#
# Daniel Drinan <ddrinan@uw.edu> 2016-12-06
#
# Script to subset loci from a GENEPOP file.
# 
# Usage:
# $> python subsetGenepop.py
#
################################################################################


import os.path
import numpy as np
import datetime, time

########################
# Start of choices
########################
choice1 = ''
while choice1 != 'A' and choice1 != 'B':
    choice1 = raw_input("Do you want to (A) select a random subset of loci or (B-untested) select a specific set of loci? ")

num_loci = 0
file_with_loci = ''
loci_to_retain = []




################################################################################
## 
## Are we retaining a random set of loci? How many?
##
################################################################################
if choice1 == 'A': # how many loci
    while num_loci < 1:
        num_loci = raw_input("How many random loci do you want to subset? ")

        try: # test to see if an integer was supplied at the command line
            num_loci = int(num_loci)
        except ValueError:
            print("That's not an integer")
            num_loci = 0


################################################################################
## 
## Are we retaining a known set of loci? Where is the file?
##
################################################################################
elif choice1 == 'B': # what is the file that has the list of loci
    while not os.path.isfile(file_with_loci):
        file_with_loci = raw_input("Which file contains the set of loci (one locus per line)? ")


    # Read in list of loci to be retained
    locus_file = open(file_with_loci, 'r')
    for line in locus_file:
        loci_to_retain.append(line.split()[0])
    locus_file.close()




################################################################################
## 
## What is the GENEPOP file to subsample from?
##
################################################################################
genepop_file = ''
while not os.path.isfile(genepop_file):
    genepop_file = raw_input("Which GENEPOP file would you like to subsample from? ")




################################################################################
## 
## How are locus names formatted in the GENEPOP file?
##
################################################################################
locus_format = ''
while locus_format != 'A' and locus_format != 'B':
    locus_format = raw_input("In the GENEPOP file, are locus names (A) comma separated on a single line or (B) separated on different lines? ")




################################################################################
## 
## Parse locus names
##
################################################################################

genepop_input = open(genepop_file)
loci = []
line = genepop_input.readline() # header line
line = genepop_input.readline() # either first locus or comma separated line of all loci

if locus_format == 'A':
    tmp_line = line.replace(' ', '')
    loci = tmp_line.split(',')

elif locus_format == 'B':
    while 'Pop' not in line and 'pop' not in line and 'POP' not in line:
        loci.append(line.split()[0])
        line = genepop_input.readline()

print 'loci =', loci


################################################################################
## 
## Check to see if there are enough loci to subsample OR that all loci to 
## subsample exist. Parse file and select the subset of loci to keep.
##
################################################################################
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
output_genepop = 'Subsample of loci ' + st + '\n'

number = []
if choice1 == 'A' and len(loci) > num_loci: # select a random set of loci
    number = np.random.choice(len(loci), num_loci, replace=False)
elif choice1 == 'B':
    tmp_counter = -1
    for locus in loci: # 'loci' is a list of loci that exist in the genepop file
        tmp_counter += 1
        if locus in loci_to_retain: # names of all the loci to keep
            number.append(tmp_counter)

for item in number:
    output_genepop += loci[item] + '\n'



################################################################################
## 
## Parse through genotypes in GENEPOP file to subsample
##
################################################################################
if choice1 == 'B': # 
    line = genepop_input.readline()

while line:
    if 'pop' in line or 'POP' in line or 'Pop' in line:
        output_genepop += line
    elif len(line) > 2:
        tmp_line = line.replace(' ,', ',').replace('\t,', ',')
        
        output_genepop += tmp_line.split(',')[0]
        output_genepop += ','
        line_split = tmp_line.split()[1:]

        for i in number:
            print i, len(line_split), tmp_line[0:30]
            output_genepop += ' ' + line_split[i] + ' '

        output_genepop += '\n'
        

    line = genepop_input.readline()

genepop_input.readline()



################################################################################
## 
## Print output file
##
################################################################################
output_file = raw_input("Where do you want to save your new GENEPOP file? ")

output_file_genepop = open(output_file, 'w')
output_file_genepop.write(output_genepop)
output_file_genepop.close()

