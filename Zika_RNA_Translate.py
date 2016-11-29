#!usr/bin/env/ python

gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
g = open("Amino_Acid.fa", "w")
f = open("FASTA.fa", "r")
content = [x.strip("\n") for x in f.readlines()]
RNA_Genome =''
for element in content:
	if element.startswith(">"):
		continue
	else:
		RNA_Genome += element


start_codon = RNA_Genome.find("ATG")
protein = ''
for nucleotide in range(start_codon, len(RNA_Genome), 3):
	codon= RNA_Genome[nucleotide:nucleotide+3]
	if len(codon) < 3:
		break	
	else:
			amino_acid = gencode[codon]
			if amino_acid == "_":
				break
			else:
				protein += amino_acid
		
g.write(protein)
g.close()




