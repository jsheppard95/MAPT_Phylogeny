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
header = {}
aligned_seqs = {}
with open(fname, "r") as f:
    line = f.readline()
    seq_header = False
    seq_data = False
    while line != "":
        line = f.readline()

        # Identify Aligned Sequence Section
        if line == "BEGIN CHARACTERS;\n":
            seq_header = True
            line = f.readline()
        elif "MATRIX" in line:
            seq_header = False
            seq_data = True
            line = f.readline()

        # Read Nexus file flags
        if seq_header:
            line_list = line.split()
            # Remove ";" from end string
            line_list[-1] = line_list[-1][:-1]

            if line_list[0] == "FORMAT":
                header[line_list[0] + line_list[1]] = line_list[3]
                header[line_list[4]] = line_list[6]
                header[line_list[7]] = line_list[9]
            else:
                header[line_list[0]] = line_list[1]

        elif seq_data:
            # Remove leading \t
            line = line[1:]
            if len(line) > 0:
                if line[0] == "'":
                    # e.g 'Homo sapiens XP_011509496'   MADER...
                    line_list = line.split("'")
                    # Remove first item - empty string
                    line_list.pop(0)
                    # Remove leading white space and \n from sequence
                    line_list[1] = line_list[1].strip()
                else:
                    # e.g Villosa_lienosa_MAP_JR498090  -----...
                    line_list = line.split()
                aligned_seqs[line_list[0]] = line_list[1]
            else:
                # finshed reading sequences
                seq_data = False

print(header)
print(aligned_seqs)