"""
Script to convert protein sequence files from FASTA (.faa, .afa) to Phylip format
"""
from Bio import AlignIO

fasta_file = "proteins_aligned.afa"
in_format = "fasta"
phy_file = "proteins_aligned.phy"
out_format = "phylip"

AlignIO.convert(fasta_file, in_format, phy_file, out_format)