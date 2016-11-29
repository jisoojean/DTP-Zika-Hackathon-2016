#!/usr/bin/python
from sys import argv
from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO

script, input_file, output_file = argv

needle_cline = NeedleCommandline(asequence="zika_ref_aaseq.fasta", bsequence=str(input_file), gapopen=10, gapextend=0.5, outfile=str(output_file))
stdout, strderr = needle_cline()

align = AlignIO.read(str(output_file), "emboss")

list_mutations = []
for i in range(align.get_alignment_length()):
	if align[0,i] == align[1,i]:
		pass
	elif align[0,i] == '-' or align[1,i] == '-':
		pass
	elif align[0,i] != align[1,i]:
		list_mutations.append(i)

print list_mutations
