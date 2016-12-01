#!/usr/bin/python
from sys import argv
from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO
from cogent import DNA
from cogent.core.genetic_code import DEFAULT as standard_code
import os

script, input_file = argv

filename = ""
for letter in input_file:
	if letter == ".":
		break
	else:
		filename += letter

f = open(input_file, 'r')
g = open("%s_aa_seq.fasta" % filename, 'w')

genome_identity = ""
genome_seq = ""
for line in f:
	if line[0] == '>':
		genome_identity = line.rstrip('\n')
	else: 
		genome_seq += line.rstrip('\n')

seq = DNA.makeSequence(genome_seq)
translated_seq = standard_code.sixframes(seq)

for frame in translated_seq:
	protein = frame.split("*")
	for peptide in protein:
		if len(peptide) > 3000:
			viral_protein = peptide

for number in range(len(viral_protein)):
	if viral_protein[number] == "M":
		viral_protein = viral_protein[number:]
		break

g.write(genome_identity + '\n')
g.write(viral_protein)
f.close()
g.close()

# Pairwise alignment by Needleman-Wunsch algorithm
needle_cline = NeedleCommandline(asequence="zika_ref_aa_seq.fasta", bsequence="%s_aa_seq.fasta" % filename, gapopen=10, gapextend=0.5, outfile="%s_pwa.txt" % filename)
stdout, strderr = needle_cline()

# Converting the alignment result into a manipulable format
align = AlignIO.read("%s_pwa.txt" % filename, "emboss")

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

os.system("mkdir %s" % filename)
os.system("mv %s_aa_seq.fasta %s_pwa.txt %s" % (filename, filename, filename))
