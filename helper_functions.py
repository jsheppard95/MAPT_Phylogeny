"""
helper_functions.py
Jackson Sheppard
"""
import requests
import json
import xmltodict
import os

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