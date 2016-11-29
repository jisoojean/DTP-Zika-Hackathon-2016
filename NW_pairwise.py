#!/usr/bin/python
from sys import argv
from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO

script, input_file = argv

needle_cline = NeedleCommandline(asequence="zika_ref_aaseq.fasta", bsequence=str(input_file), gapopen=10, gapextend=0.5, outfile="needle.txt")
stdout, strderr = needle_cline()

align = AlignIO.read("needle.txt", "emboss")

length_of_alignment = align.get_alignment_length()
list_mutations = []
list_insertions = []
list_deletions = []

for i in range(length_of_alignment):
	if align[0,i] == align[1,i]:
		pass
	elif align[0,i] == '-':
		list_insertions.append(i)
	elif align[1,i] == '-':
		list_deletions.append(i)
	elif align[0,i] != align[1,i]:
		list_mutations.append(i)
