#!/usr/bin/python
import sys
from sys import argv
import numpy
import operator
from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO
from cogent import DNA
from cogent.core.genetic_code import DEFAULT as standard_code
import matplotlib as mpl
import matplotlib.pyplot as plt
import os

sys.tracebacklimit = 0
try:
 	script, input_file = argv
except:
	raise ValueError("You need to input the FASTA file as an argument.")

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

print "Translating your input DNA sequence..."

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

print "%s_aa_seq.fasta has been created." % filename

print "Aligning the translated sequence with the reference amino acid sequence..."

# Pairwise alignment by Needleman-Wunsch algorithm
needle_cline = NeedleCommandline(asequence="zika_ref.fasta", bsequence="%s_aa_seq.fasta" % filename, gapopen=10, gapextend=0.5, outfile="%s_pwa.txt" % filename)
stdout, strderr = needle_cline()

# Converting the alignment result into a manipulable format
align = AlignIO.read("%s_pwa.txt" % filename, "emboss")

print "%s_pwa.txt has been created." % filename

# Creating lists to store information of mutations
length_of_alignment = align.get_alignment_length()
list_substitutions = []
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
		list_substitutions.append(i)

print "Generating the graph of your proteome..."
# Creating a figure and setting its dimensions and axes
fig = plt.figure(figsize=(20, 5.0))
ax1 = fig.add_axes([0.05, 0.60, 0.9, 0.15])

# Setting colour scheme of colourbar & defining no. of proteins in ZIKV
cmap = plt.get_cmap('rainbow')
protein_number = 11

# numpy.linspace(x,y,z) creates z equally spaced digits between x and y
# assign colour map value for each digit created
colors = [cmap(i) for i in numpy.linspace(0,1,protein_number)] 
# RGBA values printed if print colors: 11 colours created for 11 proteins

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

bounds = [i[0] for index,i in enumerate(sorted(protein_names.values())) if index==0]
for i in sorted(protein_names.values()):
	bounds.append(i[1])

# creating the horizontal colour bar with the defined bounds & label
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
				norm=norm,
				boundaries=bounds,
				ticks=bounds,
				spacing='proportional',
				orientation='horizontal')
cb1.set_label('Zika Virus, Encoded Proteins')

# add protein names in appropriate colour blocks in the colour bar
for protein in sorted(protein_names.items(), key=operator.itemgetter(1)):
	ax1.text(numpy.mean([protein[1][0],protein[1][1]])/float(length_of_alignment), 0.3, protein[0],
		verticalalignment='bottom', horizontalalignment='center',
		transform=ax1.transAxes,
		color='black', fontsize=12)

# For every mutation, display a * above the colorbar at the mutation site
for mutation in list_substitutions:
	ax1.text(mutation/float(length_of_alignment), 1.0, '*',
		verticalalignment='bottom', horizontalalignment='center',
		transform=ax1.transAxes,
		color='black', fontsize=10)

# For every insertion, display a + below the colorbar at the insertion site
for mutation in list_insertions:
	ax1.text(mutation/float(length_of_alignment), -0.3, '+',
		verticalalignment='bottom', horizontalalignment='center',
		transform=ax1.transAxes,
		color='black', fontsize=15)

# For every deletion, display a ^ below the colorbar at the deletion site
for mutation in list_deletions:
	ax1.text(mutation/float(length_of_alignment), -0.3, '^',
		verticalalignment='bottom', horizontalalignment='center',
		transform=ax1.transAxes,
		color='black', fontsize=15)

# Defining the contents of the legend
legend_contents = 'Legend:\n     * = substitutions\n     + = insertions\n     ^ = deletions'

# Creating the legend 
ax1.text(0, -2, legend_contents, style='normal',
	bbox={'facecolor':'white', 'alpha':1, 'pad':10})

# Save the finished plot with the legend
fig.savefig('%s_proteome.png' % filename)

print "%s_proteome.png has been created." % filename

os.system("mkdir %s" % filename)
os.system("mv %s_aa_seq.fasta %s_pwa.txt %s_proteome.png %s" % (filename, filename, filename, filename))
