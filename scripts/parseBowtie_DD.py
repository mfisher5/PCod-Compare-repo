# script written by Dan Drinan
# cd /mnt/hgfs/Shared\ Drive\ D/Pacific\ cod/DataAnalysis/UCstacksL1/bowtie
# python parseBowtie.py batch303_bowtieOut.sam batch303_bowtieOut_filtered.fa


import sys

input_file = open(sys.argv[1], 'r')

loci_to_remove = []
loci_to_keep = []
loci = {}
counter = 0

for line in input_file:
    counter += 1
    sys.stdout.write("\rNumber of Bowtie output lines read: %d" % counter)
    locus0 = line.split()[0] # specifies which columns to extract
    locus1 = line.split()[2]
    seq = line.split()[9]

    loci[locus0] = seq

    if locus0 == locus1: # keep it
        loci_to_keep.append(locus0)
    else:
        loci_to_remove.append(locus0)
        loci_to_remove.append(locus1)


input_file.close()

sys.stdout.write('\n')

# write "good" loci to a file
output_file = open(sys.argv[2], 'w')
counter = 0

for item in loci.keys():
    # in order to keep a locus, it has to have aligned only to itself.
    if (item not in loci_to_remove) and (item in loci_to_keep):
        output_file.write('>' + item + '\n')
        output_file.write(loci[item] + '\n')
        counter += 1
        sys.stdout.write("\rNumber of sequences written to output: %d" % counter)


output_file.close()
sys.stdout.write('\n')



