# Import dependencies
import requests
import os
import time
import numpy as np
import pandas as pd
import helper_functions as hf

# Read protein accession numbers from proteins.txt
infile = os.path.join("accession_files", "proteins.csv")
protein_df = pd.read_csv(infile)
print("Input Protein List:")
print(protein_df.head())

protein_df = hf.get_gids_sequences(protein_df)
print("Resulting Dataframe:")
print(protein_df.head())


# Protein Search Summary
print("Total Number of Distinct Proteins in Input File:", len(protein_df))
print("")
print("Number of Proteins From db=protein:", len(protein_df[protein_df["DB"] == "protein"]))
print("")
print ("Proteins From db=nuccore (%s):" % len(protein_df[protein_df["DB"] == "nuccore"]))
print(protein_df[protein_df["DB"] == "nuccore"])
print("")
print("Proteins Not Found (%s):" % len(hf.show_NaN_rows(protein_df)))
print(hf.show_NaN_rows(protein_df))
print("")

# Checking proteins whose sequence lengths do not match the aa counts in the input file
# Most obvious Mismatches:
# nucleotide sequences for correct protein (GL477576, BAHO01035973, KE993814, NW_003943621)
# mRNA sequence for incorrect protein (CT004140)
ret_cnts = []
for i in range(len(protein_df)):
    protein = protein_df.iloc[i]["accession_num"]
    in_cnt = protein_df.iloc[i]["aa_cnt"]  # aa count from infile
    gid = protein_df.iloc[i]["GID"]
    if gid is not np.NaN:
        fasta_seq = hf.get_fasta_from_df(protein_df, protein)
        seq = hf.convert_fasta_to_str(fasta_seq)
        ret_cnt = len(seq)
        ret_cnts.append(ret_cnt)
    else:
        ret_cnts.append(np.NaN)

print("Input vs. Returned AA Count:")
protein_df["api_returned_aa_count"] = ret_cnts

# Protein sequence counts that do not match the input count:
aa_mismatch_df = protein_df[protein_df["aa_cnt"] != protein_df["api_returned_aa_count"]]
print(aa_mismatch_df.head())  # 66 of 88 proteins have this issue, cannot drop them all!

# Building df of Proteins to be removed:
protein_NaNs = hf.show_NaN_rows(protein_df)

nuccore_df = protein_df[protein_df["DB"] == "nuccore"]

duplicate_species_df = protein_df[protein_df.duplicated("species")]  # default keep='first' -> keep first occurence

protein_drop_df = pd.concat([protein_NaNs, nuccore_df, duplicate_species_df])
print("List of Proteins To Drop:")
print(protein_drop_df)

# For now, remove the above protein sequences
# These are either:
# NaNs, i.e no returned sequence
# nucleotide sequences for correct protein (GL477576, BAHO01035973, KE993814, NW_003943621) (from nuccore)
# mRNA sequence for incorrect protein (CT004140) (from nuccore)
# Alternative sequence for the same species - currently just keeping the MAPT version

# Remove protein_drop_df (a subset of protein_df) from protein_df:
# First merge the df's with `inidcator` showing if the column is left_only, both, or right_only
protein_filt_df = hf.remove_subset_from_df(protein_df, protein_drop_df)
protein_filt_df = protein_filt_df.astype({"api_returned_aa_count": "int"})
print(protein_filt_df)
print("Number of Proteins:", len(protein_filt_df))

hf.show_NaN_rows(protein_filt_df)

# Write cleaned sequence series to file in FASTA format
protein_fasta_file = os.path.join("fasta_files", "proteins_unique.faa")
with open(protein_fasta_file, "w") as f:
    for sequence in protein_filt_df["Sequence"]:
        f.write(sequence[:-1])

# Write protein_filt_df to csv for use elsewhere
protein_csv_file = os.path.join("Protein_Sequence_Search_Output", "proteins.csv")
protein_filt_df.to_csv(protein_csv_file, index=False)

# Write nuccore sequences to csv for nucleotide translation script
protein_df[protein_df["DB"] == "nuccore"].to_csv(os.path.join("Protein_Sequence_Search_Output", "nucleotides.csv"), index=False)

# Checking MAPT proteins only
# still removing nuccore proteins (GL477576, CT004140)
# and those where sequences not found (scaffold11486, JL1528)
infile = os.path.join("accession_files", "mapt_only.txt")
with open(infile, "r") as f:
    lines = f.readlines()
mapt_only = [line.replace("\n", "") for line in lines]
print(len(mapt_only))

is_mapt = protein_df["accession_num"].isin(mapt_only)
mapt_df = protein_df[is_mapt]

mapt_fasta_file = os.path.join("fasta_files", "mapt_only.faa")
with open(mapt_fasta_file, "w") as f:
    for sequence in mapt_df["Sequence"]:
        f.write(sequence[:-1])