
infile = open("genepop/batch8/batch_8_filteredMAF_filteredIndivids30_filteredLoci_filteredHWE_filteredRepsC_INF.txt", "r")
outfile = open("batch_8_maf_by_monoVpoly.txt", "w")

loci_list = []
maf_mono_by_locus = {}
maf_poly_by_locus = {}
n_pops = 9

# skip through the INF file until summary table is reached
line = infile.readline()
while not line.startswith("Tables of allelic frequencies for each locus:"):
    line = infile.readline()

# calculate # sample sites monomorphic for each locus, using summary table
line = infile.readline() # blank line after table name
while line:
    if line.startswith(" Locus:"):
        locus = line.strip().split()[1]
        loci_list.append(locus)
        mono_at_locus_allele1 = []
        poly_at_locus_allele1 = []
        mono_at_locus_allele2 = []
        poly_at_locus_allele2 = []
        n_mono = 0
        n_poly = 0
        line = infile.readline() # extra lines
        line = infile.readline() # pop-alleles-genes header
        line = infile.readline() # more extra lines
        line = infile.readline() # alleles 1 - 2 header
        for i in range(0,n_pops):
            line = infile.readline()
            freq1 = float(line.strip().split()[1])
            freq2 = float(line.strip().split()[2])
            if freq1 == 0.0:
                mono_at_locus_allele1.append(0.0)
                mono_at_locus_allele2.append(1.0)
                n_mono += 1
            elif freq2 == 0.0:
                mono_at_locus_allele2.append(0.0)
                mono_at_locus_allele1.append(1.0)
                n_mono += 1
            else:
                poly_at_locus_allele1.append(freq1)
                poly_at_locus_allele2.append(freq2)
                n_poly += 1
        if n_mono == 0:
            maf_mono_by_locus[locus] = "NA"
        else:
            mono_a1_avg_freq = float(sum(mono_at_locus_allele1))/float(len(mono_at_locus_allele1) + 0.000000000001)
            mono_a2_avg_freq = float(sum(mono_at_locus_allele2))/float(len(mono_at_locus_allele2) + 0.000000000001)
            maf_mono_by_locus[locus] = min(mono_a1_avg_freq, mono_a2_avg_freq)
        if n_poly == 0:
            maf_poly_by_locus[locus] = "NA"
        else:
            poly_a1_avg_freq = float(sum(poly_at_locus_allele1))/float(len(poly_at_locus_allele1) + 0.000000000001)
            poly_a2_avg_freq = float(sum(poly_at_locus_allele2))/float(len(poly_at_locus_allele2) + 0.000000000001)
            maf_poly_by_locus[locus] = min(poly_a1_avg_freq, poly_a2_avg_freq)
        
    line = infile.readline()
infile.close()

print "processed ", len(loci_list), " loci."
print "writing output file..."

outfile.write("locus\tmaf_monomorhpic\tmaf_polymorphic\n")

for locus in loci_list:
    outfile.write(locus + "\t" + str(maf_mono_by_locus[locus]) + "\t" + str(maf_poly_by_locus[locus]) + "\n")

outfile.close()
print "Done."