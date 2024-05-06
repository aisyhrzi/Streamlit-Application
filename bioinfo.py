import streamlit as st
import requests
from Bio import pairwise2
from Bio.Seq import Seq
from xml.etree import ElementTree as ET
from Bio.SeqUtils import molecular_weight

# Global variable to store the protein sequence retrieved from UniProt
global_protein_sequence = ""  # Initialize as empty

# Streamlit page configuration
st.set_page_config(page_title="Protein Data Analysis", layout="wide")

def main():
    st.title("Protein Data Analysis App")
    protein_id = st.sidebar.text_input("Enter UniProt ID", value="P04637")  # Default ID for TP53 human

    analyze_button = st.sidebar.button("Analyze Protein")

    if analyze_button:
        global global_protein_sequence  # Use the global variable to store the retrieved sequence
        protein_data = fetch_protein_data(protein_id)
        if protein_data:
            global_protein_sequence = protein_data["sequence"]  # Save retrieved sequence globally
            display_protein_info(protein_data)

    st.sidebar.subheader("Protein Sequence Analysis")
    sequence_input = st.sidebar.text_area("Enter Protein Sequence", value="")
    sequence_button = st.sidebar.button("Analyze Sequence")

    if sequence_button and sequence_input and global_protein_sequence:
        analyze_protein_sequence(sequence_input)

# Function to fetch protein data and parse XML
def fetch_protein_data(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.xml"
    response = requests.get(url)
    if response.status_code == 200:
        root = ET.fromstring(response.content)
        description = root.findtext('.//{http://uniprot.org/uniprot}fullName')
        sequence = root.findtext('.//{http://uniprot.org/uniprot}sequence').replace('\n', '').strip()
        return {
            "description": description,
            "sequence": sequence
        }
    else:
        st.error("Failed to retrieve data.")
        return None

# Display protein characteristics
def display_protein_info(data):
    st.subheader("Protein Characteristics")
    st.write("Description:", data["description"])
    st.write("Protein Length:", len(data["sequence"]))
    st.write("Molecular Weight: {:.2f} Da".format(molecular_weight(data["sequence"], seq_type='protein')))

# Function to analyze a protein sequence by aligning it against the retrieved UniProt sequence
def analyze_protein_sequence(sequence):
    global global_protein_sequence  # Access the global protein sequence
    st.write("Aligning against retrieved UniProt sequence...")
    align_sequences(global_protein_sequence, sequence)

# Align two sequences and show the alignment
def align_sequences(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)  # Global alignment between the two sequences
    if alignments:
        alignment_text = pairwise2.format_alignment(*alignments[0])  # Format the first alignment result
        st.text("Alignment:")
        st.text(alignment_text)
    else:
        st.warning("No alignments found.")

if __name__ == "__main__":
    main()
