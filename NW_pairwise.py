#!/usr/bin/python
from sys import argv
from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO

script, input_file = argv

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

g = open("%s_aa_seq.fasta" % input_file , "w")
f = open(input_file, "r")

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

# Pairwise alignment by Needleman-Wunsch algorithm
needle_cline = NeedleCommandline(asequence="zika_ref_aaseq.fasta", bsequence="%s_aa_seq.fasta" % input_file, gapopen=10, gapextend=0.5, outfile="%s_pwa.txt" % input_file)
stdout, strderr = needle_cline()

# Converting the alignment result into a manipulable format
align = AlignIO.read("needle.txt", "emboss")

length_of_alignment = align.get_alignment_length()
list_mutations = []
list_insertions = []
list_deletions = []

# Finding substitution, insertions and deletions
for i in range(length_of_alignment):
	if align[0,i] == align[1,i]:
		pass
	elif align[0,i] == '-':
		list_insertions.append(i)
	elif align[1,i] == '-':
		list_deletions.append(i)
	elif align[0,i] != align[1,i]:
		list_mutations.append(i)

print list_mutations
print list_deletions
print list_insertions

# Default library of protein locations in the genome
protein_names = {
'C': (2, 104), #103
'PrM': (123, 290), #168
'E': (291, 790), #93
'NS1': (791, 1142), #352
'NS2A': (1143, 1368), #226
'NSB': (1369, 1498), #130
'NS3': (1499, 2115), #617
'NS4A': (2116, 2242), #127
'P2K': (2243, 2265), #23
'NS4B': (2266, 2516), #251
'NS5': (2517, 3419)} #903

# Modifications of the dictionary depending on the presence of insertions and deletions
for dlt in list_deletions:
	for k, v in protein_names.items():
		if dlt <= v[0]:
			protein_names[k] = (v[0]+1, v[1]+1)
		elif dlt > v[0] and dlt <= v[1]:
			protein_names[k] = (v[0], v[1]+1)

for insrt in list_insertions:
	for k, v in protein_names.items():
		if insrt <= v[0]:
			protein_names[k] = (v[0]+1, v[1]+1)
		elif insrt > v[0] and insrt <= v[1]:
			protein_names[k] = (v[0], v[1]+1)

print protein_names



