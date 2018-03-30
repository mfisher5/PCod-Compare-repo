### This file with take the output of multiple `samtools depth` files and combine into a text file ###

individs = open("Korean_subset_individs.txt", "r")
outfile = open("KOR_catalogloci_depths.txt", "w")
outfile.write("sample\tn_cat_loci\tdepths\n")


sample_list_kor = []
for line in individs:
	sample_list_kor.append(line.strip())
individs.close()

catalog_loci_dict_kor = {}
depth_dict_kor = {}
n_cat_all_kor = []
depth_all_kor = []

for sample in sample_list_kor:
	infile = open("../../../PCod-Korea-repo/stacks_b7_wgenome/" + sample + "_depth.txt", "r")
	sample_cat_loci = []
	sample_depths = []
	for line in infile:
		linelist = line.strip().split()
		if linelist[0] not in sample_cat_loci:
			sample_cat_loci.append(linelist[0])
			sample_depths.append(int(linelist[2]))
	infile.close()
	catalog_loci_dict_kor[sample] = sample_cat_loci
	depth_dict_kor[sample] = sample_depths
	outfile.write(sample + "\t" + str(len(sample_cat_loci)) + "\t" + ",".join(str(x) for x in sample_depths) + "\n")

outfile.close()
