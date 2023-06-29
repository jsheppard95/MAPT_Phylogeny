"""
phylo.py

Author: Jackson Sheppard
Last Edit: 06/19/23

Script to read aligned MAPT/MAP2/MAP4 Protein sequences for reproduction of
Bayesian consensus phylogenetic tree in Sundermann et al. BMC Genomics (2016)
17:264 Fig. 2a.
"""

from itertools import islice
import numpy as np
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

def take(n, iterable):
    """Return the first n items of the iterable as a list."""
    return list(islice(iterable, n))

print(take(10, header.items()))
#print(take(10, aligned_seqs.items()))

# Write sequences to output file in PHYLIP (*.phy) format.
print("1st Sequence Lenghth:", len(aligned_seqs["Homo sapiens XP_011509496"]))
print("2nd Sequence Lenghth:", len(aligned_seqs["Nomascus leucogenys XP_012357075"]))
print("10th Sequence Lenghth:", len(aligned_seqs["Orcinus orca XP_004262822"]))
seq_length = len(aligned_seqs["Homo sapiens XP_011509496"])
# Clean species names:
# Consistent name format, spaces replaced by underscores
species = list(aligned_seqs.keys())
for id in species:
    if " " in id:
        new_id = id.replace(" ", "-")
        aligned_seqs[new_id] = aligned_seqs.pop(id)
    elif "_" in id:
        new_id = id.replace("_", "-")
        aligned_seqs[new_id] = aligned_seqs.pop(id)
#print(take(30, aligned_seqs.items()))

# Check for unique species names for PHYLIP file
SPECIES_NAME_LEN = 10
species = list(aligned_seqs.keys())
phyl_ids = []
for id in species:
    id_spl = id.split("-")
    phyl_ids.append(id_spl[-1][:SPECIES_NAME_LEN])
phyl_ids = np.array(phyl_ids)
print(len(phyl_ids))
print(len(np.unique(phyl_ids)))
unique_ids = len(phyl_ids) == len(np.unique(phyl_ids))
print("Unique:", unique_ids)

# Confirmed unique keys, update sequence dictionary
for id in species:
#    print(id)
    id_spl = id.split("-")
    new_id = id_spl[-1][:SPECIES_NAME_LEN]
    aligned_seqs[new_id] = aligned_seqs.pop(id)

# Write aligned sequences in PHYLIP format - Using Sequential
# 5 42 
# Turkey                  AAGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGGTAT 
# Salmo_schiefermuelleri  AAGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGTAT 
# H_sapiens               ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA 
# Chimp                   AAACCCTTGCCGTTACGCTTAAACCGAGGCCGGGACACTCAT 
# Gorilla                 AAACCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCTTAA
# 5 = # of species, 42 = length of aligned sequences
# species has 10 character limit - use first 10 characters (unique)

outfile = os.path.join("TreeBASE", "MAPPHY.phy")
with open(outfile, "w") as f:
    f.write(str(len(aligned_seqs)))
    f.write(" " + str(seq_length))
    f.write(" \n")
    for id in aligned_seqs.keys():
#        print(id)
        if len(id) < SPECIES_NAME_LEN:
            f.write(id + (SPECIES_NAME_LEN-len(id))*" ")
        else:
            f.write(id)
        f.write(aligned_seqs[id] + "\n")
