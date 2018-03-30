# Code by Eleni, November 28, 2016
# unix code to match a string in a directory of files. I wrote this because I wanted to count the number of loci present in each individual sample after pstacks. 
# My assumption here is that the word "consensus" appears once for each reported locus in the individual.tags.tsv files produced by pstacks. So if I count the 
# occuences of that word, I count the number of loci present for each individual sample. 


#SYNOPSIS
#       grep [OPTIONS] PATTERN [FILE...]

grep --count --with-filename consensus *.tags.tsv > output_file.txt

#	explanation of terms:
#		grep = print lines matching a pattern
#		--count =  Suppress normal output; instead print a count of matching  lines  for  each  input  file.
#		--with-filename  = Print  the  file  name for each match.  This is the default when  there is more than one file to search.
#		consensus = PATTERN to be matched ("consensus" is the pattern in this example)
#		*.tags.tsv = FILE; do the grep stuff in any file in the current directory that ends with the extension .tags.tsv
#		> output_file.txt = append the results of code to this output file
