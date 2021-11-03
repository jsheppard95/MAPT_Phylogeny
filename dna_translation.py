# Script to read in DNA sequence and convert to amino acid sequence
import os
import pandas as pd
import requests
from helper_functions import api_call, get_fasta_from_df, convert_fasta_to_str

# Read protein accession numbers from proteins.txt
infile = os.path.join("Protein_Sequence_Search_Output", "nucleotides.csv")
nucleotide_df = pd.read_csv(infile)

# Adding one more nucleotide for testing:
mus = "NM_207618"
mus_data = api_call("nuccore", mus)
gid = mus_data["eSummaryResult"]["DocSum"]["Id"]
gid_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sequences&id=" + gid + "&rettype=fasta&retmode=text"
gid_resp = requests.get(gid_url)
seq = gid_resp.content.decode("utf-8")
seq
nucleotide_df = nucleotide_df.append({
    "accession_num" : mus,
    "GID": gid,
    "Sequence": seq
}, ignore_index=True)

fasta_seq = get_fasta_from_df(nucleotide_df, mus)
seq = convert_fasta_to_str(fasta_seq)
ret_cnt = len(seq)
nucleotide_df.loc[5, "api_returned_aa_count"] = ret_cnt
print(nucleotide_df)
print("")

# Each amino acid created from 3 nucleotides (triplet)
# according to the mapping below
def translate(dna_seq):
    """
    Function to convert DNA nucleotide sequences into
    amino acid protein sequences
    
    Parameters:
    -----------
    dna_seq : str
        DNA sequence as string, e.g: TATATATATGTAAGGTT...
    """
    translation_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    aa_seq = ""
    # Check DNA sequence is divisible by 3, otherwise throw an error
    if len(dna_seq) % 3 == 0:
        # Massage NCBI sequence to pure string of nucleotides,
        # TATATATATGTA...
        for i in range(0, len(dna_seq), 3):
            codon = dna_seq[i:i + 3]
            aa_seq += translation_table[codon]
        return aa_seq
    else:
        print("Invalid DNA sequence, length not divisible by 3.")


CDS_START = 21  # Coding region start, NB: nucleotides numbered from 1
CDS_END = 938  # Coding region end
print("DNA Sequence:", mus)
print(seq)
print("")
print("Translated Protein Sequence:")
translated_seq = translate(seq[CDS_START - 1: CDS_END])
print(translated_seq)
print("")
actual_translation = "MSTHDTSLKTTEEVAFQIILLCQFGVGTFANVFLFVYNFSPISTGSKQRPRQVILRHMAVANALTLFLTIFPNNMMTFAPIIPQTDLKCKLEFFTRLVARSTNLCSTCVLSIHQFVTLVPVNSGKGILRASVTNMASYSCYSCWFFSVLNNIYIPIKVTGPQLTDNNNNSKSKLFCSTSDFSVGIVFLRFAHDATFMSIMVWTSVSMVLLLHRHCQRMQYIFTLNQDPRGQAETTATHTILMLVVTFVGFYLLSLICIIFYTYFIYSHHSLRHCNDILVSGFPTISPLLLTFRDPKGPCSVFFNC"
print("Actual Protein Sequence:")
print(actual_translation)