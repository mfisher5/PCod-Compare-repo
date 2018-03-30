### This is part two of the script that will convert a 2x2 matrix file into a genepop file. It is automatically generated with part 1. ASSUMPTION: SAMPLES ARE COLUMN HEADERS #####


infile = open('../../stacks_b8_wgenome_r05/batch_8_filteredMAF_filteredIndivids30_filteredLoci30_filteredHWE.txt', 'r')
genepop = open('../../stacks_b8_wgenome_r05/batch_8_filteredMAF_filteredIndivids30_filteredLoci30_filteredHWE_genepop.txt', 'w')
genepop.write('Korea + Alaska Data Final Genepop, stacks batch 8 final, uncorrected m3 filtered (MAF / Loci / Individs / HWE)\r\n')
print 'transposing genotypes matrix...'
data_matrix = []
for line in infile:
	tmp_line = ''
	tmp_line += line
	data_matrix.append(tmp_line.split('\t'))
infile.close()

transposed = zip(*data_matrix)

print 'writing loci into genepop file...'
locilist = transposed[0]
LociIndex = range(0, len(locilist))
for i in LociIndex:
	if transposed[0][i] != 'sample':
		genepop.write(transposed[0][i] + '\r\n')

Pohang15 = transposed[1:35]
Geoje15 = transposed[35:72]
Namhae15 = transposed[72:91]
YellowSea16 = transposed[91:121]
Jukbyeon07 = transposed[121:158]
JinhaeBay07 = transposed[158:203]
JinhaeBay08 = transposed[203:250]
Boryeong07 = transposed[250:274]
Geoje14 = transposed[274:309]
Kodiak03 = transposed[309:355]
Adak06 = transposed[355:398]
WashCoast05 = transposed[398:442]
HecStrait04 = transposed[442:489]
SalishSea12_13 = transposed[489:531]
JuandeFuca12 = transposed[531:553]
PWSound12 = transposed[553:600]
UnimakPass03 = transposed[600:644]
last_line = list(transposed[644])
seq = range(0, len(last_line))
for i in seq:
	last_line[i] = last_line[i].strip('\r\n')


print 'writing genotypes into genepop file by population...'
genepop.write('Pop' + '\r\n')
for line in Pohang15:
	linestr = '\t'.join(line[1:])
	newline = line[0] + ',\t' + linestr + '\r\n'
	genepop.write(newline)

genepop.write('Pop' + '\r\n')
for line in Geoje15:
	linestr = '\t'.join(line[1:])
	newline = line[0] + ',\t' + linestr + '\r\n'
	genepop.write(newline)

genepop.write('Pop' + '\r\n')
for line in Namhae15:
	linestr = '\t'.join(line[1:])
	newline = line[0] + ',\t' + linestr + '\r\n'
	genepop.write(newline)

genepop.write('Pop' + '\r\n')
for line in YellowSea16:
	linestr = '\t'.join(line[1:])
	newline = line[0] + ',\t' + linestr + '\r\n'
	genepop.write(newline)

genepop.write('Pop' + '\r\n')
for line in Jukbyeon07:
	linestr = '\t'.join(line[1:])
	newline = line[0] + ',\t' + linestr + '\r\n'
	genepop.write(newline)

genepop.write('Pop' + '\r\n')
for line in JinhaeBay07:
	linestr = '\t'.join(line[1:])
	newline = line[0] + ',\t' + linestr + '\r\n'
	genepop.write(newline)

genepop.write('Pop' + '\r\n')
for line in JinhaeBay08:
	linestr = '\t'.join(line[1:])
	newline = line[0] + ',\t' + linestr + '\r\n'
	genepop.write(newline)

genepop.write('Pop' + '\r\n')
for line in Boryeong07:
	linestr = '\t'.join(line[1:])
	newline = line[0] + ',\t' + linestr + '\r\n'
	genepop.write(newline)

genepop.write('Pop' + '\r\n')
for line in Geoje14:
	linestr = '\t'.join(line[1:])
	newline = line[0] + ',\t' + linestr + '\r\n'
	genepop.write(newline)

genepop.write('Pop' + '\r\n')
for line in Kodiak03:
	linestr = '\t'.join(line[1:])
	newline = line[0] + ',\t' + linestr + '\r\n'
	genepop.write(newline)

genepop.write('Pop' + '\r\n')
for line in Adak06:
	linestr = '\t'.join(line[1:])
	newline = line[0] + ',\t' + linestr + '\r\n'
	genepop.write(newline)

genepop.write('Pop' + '\r\n')
for line in WashCoast05:
	linestr = '\t'.join(line[1:])
	newline = line[0] + ',\t' + linestr + '\r\n'
	genepop.write(newline)

genepop.write('Pop' + '\r\n')
for line in HecStrait04:
	linestr = '\t'.join(line[1:])
	newline = line[0] + ',\t' + linestr + '\r\n'
	genepop.write(newline)

genepop.write('Pop' + '\r\n')
for line in SalishSea12_13:
	linestr = '\t'.join(line[1:])
	newline = line[0] + ',\t' + linestr + '\r\n'
	genepop.write(newline)

genepop.write('Pop' + '\r\n')
for line in JuandeFuca12:
	linestr = '\t'.join(line[1:])
	newline = line[0] + ',\t' + linestr + '\r\n'
	genepop.write(newline)

genepop.write('Pop' + '\r\n')
for line in PWSound12:
	linestr = '\t'.join(line[1:])
	newline = line[0] + ',\t' + linestr + '\r\n'
	genepop.write(newline)

genepop.write('Pop' + '\r\n')
for line in UnimakPass03:
	linestr = '\t'.join(line[1:])
	newline = line[0] + ',\t' + linestr + '\r\n'
	genepop.write(newline)

linestr = '	'.join(last_line[1:])
newline = last_line[0] + ',\t' + linestr + '\r\n'
genepop.write(newline)

genepop.close()

print 'done.'