# Import dependencies
import pandas as pd
import os
from Bio import SeqIO

# Read in LIB Sequence file
infile = os.path.join("MAPT_Morgan", "Morgan_MAP-phylo_Proteins-296.LIB")
record_dict = SeqIO.to_dict(SeqIO.parse(infile, "fasta"))
print("Number of Records:", len(record_dict))
print("MAPTHU1 Sequence:")
print(record_dict["MAPTPHU1"].seq)
print("MAPTHU1 Desc:")
print(record_dict["MAPTPHU1"].description)
print("MAPTHU1 Length:", len(record_dict["MAPTPHU1"].seq))
print("")

# All available attributes of `Bio.SeqRecord`
print("Example Records for MAP2PTA1")
print("id:", record_dict["MAP2PTA1"].id)
print("seq:", record_dict["MAP2PTA1"].seq)
print("name:", record_dict["MAP2PTA1"].name)
print("description:", record_dict["MAP2PTA1"].description)
print("dbxrefs:", record_dict["MAP2PTA1"].dbxrefs)
print("features:", record_dict["MAP2PTA1"].features)
print("annotations:", record_dict["MAP2PTA1"].annotations)
print("letter_annotations", record_dict["MAP2PTA1"].letter_annotations)
print("")

# Add LIB data to df
ids = []
names = []
seqs = []
descs = []

for record in record_dict:
    ids.append(record_dict[record].id)
    names.append(record_dict[record].name)
    seqs.append(record_dict[record].seq)
    descs.append(record_dict[record].description)

mapt_morgan_df = pd.DataFrame({
    "ID": ids,
    "Name": names,
    "Sequence": seqs,
    "Description": descs
})
print("LIB as pandas df:")
print(mapt_morgan_df.head())
print("")
# Accessing a sequence
print("Sequence of First Entry (MAP2HSA4)")
print(mapt_morgan_df["Sequence"][0])
print("")

# Many more records here (296) than were in paper (88 listed accession numbers)
# Which proteins in mapt_morgan_df have some correspondence to those in the paper?

# Load proteins.csv: accession_num,aa_cnt,species from the paper
protein_acc_csv = os.path.join("accession_files", "proteins.csv")
protein_paper_df = pd.read_csv(protein_acc_csv)
print("Proteins Included in Paper:")
print(protein_paper_df.head())
print("")

# Look in the Description column of mapt_morgan_df - some accession numbers listed here
print("Accesion Numbers from LIB Description Column")
print(mapt_morgan_df["Description"].head())
print("")

print("Extracting and Appending Accession Numbers From LIB File:")
# Get dataframe of proteins listing an accession number, append accession number
ids = []
names = []
seqs = []
descs = []
acc_nums = []
for index, row in mapt_morgan_df.iterrows():
    if "|ref|" in row["Description"]:
        ids.append(row["ID"])
        names.append(row["Name"])
        seqs.append(row["Sequence"])
        descs.append(row["Description"])
        acc_num = row["Description"].split("|")[3]
        acc_nums.append(acc_num)

mapt_morgan_acc_df = pd.DataFrame({
    "ID": ids,
    "Name": names,
    "Sequence": seqs,
    "Description": descs,
    "Accession_Number": acc_nums
})

# Clean up Accession_Number column
for index, row in mapt_morgan_acc_df.iterrows():
    row["Accession_Number"] = row["Accession_Number"].split(" ")[0]
acc_nums = list(mapt_morgan_acc_df["Accession_Number"])
acc_num_bases = [acc.split(".")[0] for acc in acc_nums]
mapt_morgan_acc_df["Accession_Number_Base"] = acc_num_bases
print(mapt_morgan_acc_df.head())

# Which accession numbers from the .LIB file made it into the paper?
pd.merge(mapt_morgan_acc_df, protein_paper_df, how="inner", left_on="Accession_Number_Base", right_on="accession_num")
