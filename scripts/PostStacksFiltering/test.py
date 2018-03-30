#--- count missing data in each population


PopList = ["Pohang", "Geoje"]

for pop in PopList: 
	newstr = "#next pop\n\t\tCount_MissingGenotypesByLocus_" + pop + " = float(" + pop + ".count('0000'))" + "\n\t\tNumberOf_"
	print newstr
#	newpop = str(pop)
#	first = "#next pop" + "\n\t\tCount_MissingGenotypesByLocus_" + newpop + " = float(" + newpop + ".count("0000"))" + "\n\t\tNumberOf_" + newpop + "_individuals = float(len(" #+ newpop + "))\n" 
#	second = "\t\tPercent_MissingData_" + newpop + " = float(Count_MissingGenotypesByLocus_" + newpop + "/NumberOf_" + newpop + "_individuals)\n\t\t" + "overall_percent_missingdata.append(Percent_MissingData_" + newpop + ")\n"
#	newstr = first + second
#	print newstr

