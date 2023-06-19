"""
phylo.py

Author: Jackson Sheppard
Last Edit: 06/19/23

Script to read aligned MAPT/MAP2/MAP4 Protein sequences for reproduction of
Bayesian consensus phylogenetic tree in Sundermann et al. BMC Genomics (2016)
17:264 Fig. 2a.
"""

import os

fname = os.path.join("TreeBASE", "1_1458651913_MAP-102x1953_ExaBayesConsensusEMR.nexorg")

# Read Nexus file and extract aligned sequences:
# Want this part of the file:
# Begin CHARACTERS;
#     TITLE  Microtubule_associated_protein_Tau;
#     DIMENSIONS  NCHAR=1953;
#     FORMAT DATATYPE = Protein GAP = - MISSING = ?;
#     MATRIX
#     <Species> <Accession Number>    <Aligned Sequence>
#     ...
#
# ;
with open(fname, "r") as f:
    line = f.readline()
    seq_header = False
    seq_data = False
    while line != "":
        line = f.readline()
        if line == "BEGIN CHARACTERS;\n":
            seq_header = True
        elif "MATRIX" in line:
            seq_header = False
            seq_data = True
        if seq_header:
            line_list = line.split()
            print(line_list)
        elif seq_data:
            print(line[:79])
        if line == "END;\n":
            seq_data = False

