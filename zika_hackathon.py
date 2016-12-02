#!/usr/bin/python
import sys
from sys import argv
import os
import urllib2
import time
import re
import operator
from cogent.parse.fasta import MinimalFastaParser
from cogent import LoadSeqs,DNA
from cogent.core.genetic_code import DEFAULT as standard_code
from cogent.app.muscle import align_unaligned_seqs
from cogent.app.fasttree import build_tree_from_alignment
from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO
import numpy
import matplotlib as mpl
mpl.use('Qt4Agg')
import matplotlib.pyplot as plt
from ete3 import Tree, TreeStyle, NodeStyle

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

newseqid = str(raw_input ('Please input the new sequence ID: '))
newseqcountry = str(raw_input ('Please input the country and region of origin: '))
newseqyear = str(raw_input ('Which year was it sequenced in? '))
newinfostr = ">"+newseqid+"|"+newseqcountry+"|"+newseqyear+'\n'

accessionid = []

start_time = time.time()
ncbiupdate = raw_input ('Do you want to extract all ZIKV sequences from the NCBI database?(Y/N) ').upper()
cycle = 0

if ncbiupdate == 'Y':
	os.system("rm zikv_database.fa")
	fastafile = open('zikv_database.fa','a+')
	sourcecode = urllib2.urlopen('https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/Database/nph-select.cgi?cmd=show_query&country=any&genregion=any&go=database&host=Human&isolation=any&searchin=sequence&sequence=N&sonly=on&srcfilter_labs=include&taxid=64320')
	content = sourcecode.read()
	accession_list = re.findall('[A-Z][A-Z][0-9][0-9][0-9][0-9][0-9][0-9]', content)
	print ('Found %s ZIKV complete sequenced genomes in the NCBI database.' % (len(set(accession_list))))
	print ('Downloading the Fasta files..')
	for accession in accession_list:
		if accession not in accessionid:
			accessionid.append(accession)
			genbanksource = urllib2.urlopen('https://www.ncbi.nlm.nih.gov/nuccore/%s.1?report=gilist&log$=seqview&format=text' % (accession))
			genbankcontent = genbanksource.read()
			idcode = re.findall('[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]?', genbankcontent)
			idcodestr = ''.join(map(str, idcode))
			gbsource = urllib2.urlopen('https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=genbank&sort=&id=%s&from=begin&to=end&extrafeat=976&maxplex=3' % (idcodestr))
			gbcontent = gbsource.read()
			gbfile = open('fasta.gbk', 'w+')
			gbfile.write(gbcontent)
			gbfile.seek(0)
			for line in gbfile:
				for part in line.split():
					if "country=" in part:
						x = line.split('"')
						name = x[1]
						nocolon=name.split(":")[0]
						countryname=nocolon.replace(" ","_")
					if "collection_date=" in part:
						y = line.split('"')
						collectiondate = y[1].split('-')
						for item in collectiondate:
							if len(item) == 4 and item.isdigit() == True:
								collectionyear = item
			seqident = ">"+str(accessionid[-1])+"|"+countryname+"|"+collectionyear+'\n'
			fastafile.write(seqident)
			fastasource = urllib2.urlopen('https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&sort=&id=%s&from=begin&to=end&extrafeat=976' % (idcodestr))
			fastacontent = fastasource.readlines()
			for line in fastacontent:
				if line[0]!=">":
					fastafile.write(line.rstrip())
			fastafile.write("\n")
			gbfile.truncate()
			cycle += 1
			if cycle == 3:
				break
			print ('Fasta file ')+(accessionid[-1])+(' downloaded; ')+str((100*(cycle/(float(len(set(accession_list)))))))+'%'+' completed..'

	finaltime = round((((((time.time() - start_time)))/60)), 2)
	print ('Downloaded all %s fasta files. ' % str(len(set(accession_list))))+('Total time: '+str(finaltime)+' minutes.')
	os.remove('fasta.gbk')
	gbfile.close()
	fastafile.close()

fastafile = open('zikv_database.fa','a+')

fastafile.write(newinfostr)
inputfile = open(input_file, 'r+')
for line in inputfile.readlines():
    if line[0]!=">":
        fastafile.write(line.rstrip())
fastafile.write("\n")
fastafile.close()

FastaFile = open("zikv_database.fa", 'r')
InputSeq = MinimalFastaParser(FastaFile)
Seq = LoadSeqs(data=InputSeq,moltype=DNA,aligned=False, label_to_name=lambda x: x.split("|")[0]+x.split("|")[1][0:10]+x.split("|")[2][0:4])

FastaFile.close()

print ('The input sequence has been appended to the database FASTA file.')

print ('Aligning all sequences..')

AlignedSeq = align_unaligned_seqs(Seq,DNA)
print ('Building tree..')
tree = build_tree_from_alignment(AlignedSeq,DNA)

e = open("%s_phylotree.tre" % filename, 'w')

e.write(tree.getNewick(with_distances=True))

e.close()

print "%s_phylotree.tre has been created." % filename

t = Tree("%s_phylotree.tre" % filename)

nodename = newseqid+newseqcountry[0:10]+newseqyear[0:4]

n = t.search_nodes(name=nodename)

midpoint = t.get_midpoint_outgroup()

t.set_outgroup(midpoint)

ts = TreeStyle()
ts.show_leaf_name = True
ts.mode = "c"
ts.arc_start = -180
ts.arc_span = 360

nst = NodeStyle()
nst["bgcolor"] = "LightSteelBlue"

n[0].set_style(nst)

t.render("%s_phylotree.png" % filename, w = 500, units="mm", tree_style=ts)

FastaFile.close()
e.close()

f = open(input_file, 'r')
g = open("%s_aa_seq.fasta" % filename, 'w')

genome_identity = ""
genome_seq = ""
for line in f:
	if line[0] == '>':
		genome_identity = line.rstrip('\n')
	else:
		genome_seq += line.rstrip('\n')

print "Translating your input DNA sequence.."

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

print "Aligning the translated sequence with the reference amino acid sequence.."

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

print "Generating the graph of your proteome.."

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
os.system("mv %s_aa_seq.fasta %s_pwa.txt %s_proteome.png %s_phylotree.tre %s_phylotree.png %s" % (filename, filename, filename, filename, filename, filename))
