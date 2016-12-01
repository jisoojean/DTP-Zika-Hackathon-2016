# Real-time genetic analysis of Zika virus sequences

## Introduction
Primarily transmitted during bloodmeal by Aedes aegypti mosquitoes vector, Zika virus (ZIKV) belongs to the genus Flavivirus and is classified by sequence analysis into two main genotypic groups, African and Asian. 

Since 2015, 69 countries have reported cases of Zika virus infection, with the majority in the Americas, but also parts of Asia. Most shockingly, ZIKV infection of pregnant women in the recent outbreak has been linked to head and brain malformations, including microcephaly in children at birth.

In the context of the ongoing ZIKV epidemic, we have developed a program to analyse a given ZIKV genome sequence (e.g. from samples collected during the outbreak) with respect to the existing ZIKV genome sequences from the NCBI database, providing an insight into the specific features of the strain of interest. Furthermore, it can be considered a valuable tool as it contributes to better understanding of the origin and evolution of ZIKV and it provides a framework for analysing future viral outbreaks. 


## Aim
This document shows how to use the program (**Nameofpythoncode.py**) to analyse a ZIKV genome sequence using database interrogation and phylogenetic tools, as well as genetic and proteomics analysis.

The aims of the program are:
- To compare an input ZIKV genome sequence with the existing ZIKV genome sequences available at the [NCBI Virus Variation resource](https://www.ncbi.nlm.nih.gov/genome/viruses/variation/Zika/) (*Brister et al, 2014*) in real time by multiple sequence alignment (MSA) using MUSCLE (*Edgar and Robert, 2004*).
- To subsequently perform a maximum likelihood phylogenetic analysis.
- To obtain the amino acid sequence from the input ZIKV genome sequence.
- To compare the generated amino acid sequence with respect to a [reference ZIKV coding sequence]( https://www.ncbi.nlm.nih.gov/nuccore/226377833?report=graph&tracks=[key:sequence_track,name:Sequence,display_name:Sequence,id:STD1,category:Sequence,annots:Sequence,ShowLabel:false,shown:true,order:1][key:gene_model_track,name:Genes,display_name:Genes,id:STD3,category:Genes,annots:Unnamed,Options:ShowAll,SNPs:true,CDSProductFeats:true,ShowLabelsForAllFeatures:true,NtRuler:true,AaRuler:false,HighlightMode:2,shown:true,order:3][key:feature_track,name:Other%20features---3%27UTR,display_name:3%27UTR%20Features,id:STD4,subkey:3%27UTR,category:Features,subcategory:3%27UTR%20Features,annots:Unnamed,Layout:Adaptive,LinkedFeat:Packed,shown:true,order:4][key:feature_track,name:Other%20features---5%27UTR,display_name:5%27UTR%20Features,id:STD5,subkey:5%27UTR,category:Features,subcategory:5%27UTR%20Features,annots:Unnamed,Layout:Adaptive,LinkedFeat:Packed,shown:true,order:5]&appname=ncbientreznuccore&assm_context=GCF_000882815.1&color=0&label=0&decor=0&spacing=0&v=1:10794&c=C0C0C0&gflip=false&select=null&slim=0) by pairwise alignment using the Needleman-Wunsch algorithm.
- To map mutations in the input ZIKV genome sequence by creating a graphical representation with the protein-coding sections from the reference genome.

## Prerequisities
The following software modules are required to run the code (**reorder according to code**):
- Python version 2.7
- Numerical Python (NumPy)
- Biological Python (BioPython) 
- Python Plotting (MatPlotLib)
Besides, both the ZIKV genome sequence of interest that is to be analysed and the reference ZIKV coding sequence (*zika_ref.fasta*) should be located in the working directory in FASTA format (*.fa* or *.fasta*).


## Running
Defaults in the code are setup so that if the following is run in a terminal:
`Nameofpythoncode.py` it will prompt the user to specify the name of the file of the ZIKV genome sequence of interest as an argument.
Then, the user will be asked whether the program should update the genome sequences in the database available at the [NCBI Virus Variation resource](https://www.ncbi.nlm.nih.gov/genome/viruses/variation/Zika/). If the updating option is selected, all ZIKV genome sequences will start to be downloaded. Running time for this step will be of around 4 minutes on a modern laptop.
Total running time for the code should be of approximately **x minutes** on a modern laptop.

## Expected output
Once the program is run, a new empty directory will be created in the working directory, called *filename*. This directory will contain all the output files that have been generated, as follows:
- Complete set of existing ZIKV genome sequences available at the NCBI Virus Variation resource: *zikv_database.fa*
- Interactive phylogenetic tree in the ETE Tree Browser. 
- Graphical representation of the phylogenetic tree containing the input ZIKV genome sequence and the existing ZIKV genome sequences (*zikv_database.fa*): *filename_phylotree.png*
- Amino acid sequence generated as of the input ZIKV genome sequence: *filename_aa_seq.fasta*
- Pairwise alignment of the input ZIKV amino acid sequence (*filename_aa_seq.fasta*) and the reference coding sequence (*zika_ref.fasta*): *filename_pwa.txt*
- Graphical representation of the mutations obtained in the pairwise alignment (*filename_pwa.txt*): *filename_ZIKV.png*


## License
Please see the document called LICENSE.md.


## References
1. Brister JR, Bao Y, Zhdanov SA, Ostapchuck Y, Chetvernin V, Kiryutin B, Zaslavsky L, Kimelman M, Tatusova TA (2014). Virus Variation Resource - recent updates and future directions. Nucleic Acids Res. 42(Database issue):D660-5.
2. Edgar, Robert C. (2004), MUSCLE: multiple sequence alignment with high accuracy and high throughput, Nucleic Acids Research 32(5), 1792-97.

