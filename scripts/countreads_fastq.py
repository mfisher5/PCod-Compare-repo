### This script counts the read depth of each process_radtags fastq output file ###


import argparse 

parser = argparse.ArgumentParser(description="count number of consensus seqs in .tags files")

parser.add_argument("-s", "--samples", help="file with list of samples (if paired end, should include .1 at end)")
parser.add_argument("-d", "--directory", help="local path with directory with fastq files")
parser.add_argument("-os", "--outputshell", help="bash shell script file name. must have file extension .sh")
parser.add_argument("-of", "--output", help="text file name to store read depths, with local path. Do not include file extension")

args = parser.parse_args()
import subprocess

sample_list = []
reads_list = []

samplefile= open(args.samples, "r")
outfile = open(args.outputshell, "w")
outfile.write("#!/bin/bash\n\n")
outfile.write("cd " + args.directory + "\n\n")

for line in samplefile:
	sample = line.strip().split()[0]
	sample_list.append(sample)
	newstr = 'echo ' + sample + '\nCOUNT="$(cat ' + sample + '.fq | wc -l)"\necho $(( $COUNT / 4 )) >> ../' + args.output + '_temp.txt\n\n'
	outfile.write(newstr)

samplefile.close()
outfile.close()

print "Bash shell script created"

shellfile = "./" + args.outputshell

subprocess.call(['chmod', '+x', args.outputshell])
print "Calling shell script..."
subprocess.call([shellfile])


readsfile = open("../" + args.output + "_temp.txt", "r")
finalout = open("../" + args.output + ".txt", "w")

for line in readsfile: 
	reads_list.append(line.strip())

readsfile.close()

if len(reads_list) == len(sample_list):
	for i in range(0, len(reads_list)):
		finalout.write("\n" + sample_list[i] + "\t" + reads_list[i])

finalout.close()

