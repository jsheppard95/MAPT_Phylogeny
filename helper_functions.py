"""
helper_functions.py
Jackson Sheppard
"""
import requests
import json
import xmltodict
import os
import time
import numpy as np

# Function to perform API call and return Python Dictionary containing data
def api_call(db, protein):
    """
    Performs API call to NCBI for specified database and protein
    Parameters:
    -----------
    db : str
        NCBI database to search for, i.e. 'protein', 'nuccore', etc.
    protein : str
        Protein accession number to search for
    """
    # Create url for API call
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=" + db + "&id=" + protein
    
    # Define tempory XML and JSON files
    xml_file = "data.xml"
    json_file = "data.json"

    # Perform API call
    resp = requests.get(url)
    
    # Save XML result to temporary file
    with open(xml_file, "wb") as f:
        f.write(resp.content)

    # Read XML file and convert to OrderedDict using xmltodict
    with open(xml_file, "r") as f:
        data_dict = xmltodict.parse(f.read())
    
    # Write OrderedDict to JSON file
    json_data = json.dumps(data_dict)
    with open(json_file, "w") as f:
        f.write(json_data)
    
    # Read in JSON file into regular Python dictionary
    with open(json_file, "r") as f:
        data = json.load(f)
    
    # Clean up temporary files
    os.remove(xml_file)
    os.remove(json_file)

    return data


def get_gids_sequences(protein_df):
    """
    Function to get GID numbers and Sequences for protein accession numbers
    contained in protein_df and append this data as new colums to the same
    dataframe
    """
    gids = []
    dbs = []
    seqs = []
    for i in range(len(protein_df)):
        print("Protein %s of %s" % (i + 1, len(protein_df)), end="\r")
        # Wait 1 sec every 3 calls to not bog down servers
        if (i+1) % 3 == 0:
            time.sleep(1)
        protein = protein_df.iloc[i]["accession_num"]
        # First search NCBI DB=protein
        db = "protein"
        result = api_call(db, protein)
        try:
            gid = result["eSummaryResult"]["DocSum"]["Id"]
        except KeyError:
            # if error, search NCBI DB=nuccore
            try:
                db = "nuccore"
                result = api_call(db, protein)
                gid = result["eSummaryResult"]["DocSum"]["Id"]
            # If still error, return None
            except KeyError:
                db = None
                gid = np.NaN
        if gid is not np.NaN:
            gid_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sequences&id=" + gid + "&rettype=fasta&retmode=text"
            seq = requests.get(gid_url).content.decode("utf-8")
            seqs.append(seq)
        else:
            seqs.append(np.NaN)
        gids.append(gid)
        dbs.append(db)
    protein_df["DB"] = dbs
    protein_df["GID"] = gids
    protein_df["Sequence"] = seqs
    return protein_df


def show_NaN_rows(df):
    is_NaN = df.isnull()
    row_has_NaN = is_NaN.any(axis=1)
    df_NaN = df[row_has_NaN]
    return df_NaN


# Helper function to get a protein sequence and its length
def get_sequence_and_count(protein_df, protein):
    protein_indeces = protein_df.index
    idx_array = protein_indeces[protein_df["accession_num"] == protein]
    idx = idx_array[0]
    seq_raw = protein_df[protein_df["accession_num"] == protein]["Sequence"][idx]
    seq_array = seq_raw.split("\n")
    seq = ""
    for i in range(1, len(seq_array)):
        seq += seq_array[i]
    return (len(seq), seq_raw, seq)