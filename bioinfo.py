import streamlit as st
import requests
import networkx as nx
import matplotlib.pyplot as plt
from Bio.SeqUtils import molecular_weight
from Bio import pairwise2
from Bio.Seq import Seq
from xml.etree import ElementTree as ET
import time

# Global variable to store the protein sequence fetched from UniProt in FASTA format
global_fasta_sequence = ""

# Streamlit Page Config
st.set_page_config(page_title="Protein Data Analysis", layout="wide")


def main():
    st.title("Protein Data Analysis App")
    protein_id = st.sidebar.text_input("Enter UniProt ID", value="P04637")  # Default ID for TP53 human

    analyze_button = st.sidebar.button("Analyze Protein")

    if analyze_button:
        show_progress_bar()
        protein_data = fetch_protein_data(protein_id)
        if protein_data:
            global global_fasta_sequence
            global_fasta_sequence = fetch_fasta(protein_id)
            display_protein_info(protein_data)
            display_ppi_network(protein_id)

    st.sidebar.subheader("Protein Sequence Analysis")
    sequence_input = st.sidebar.text_area("Enter Protein Sequence", value="")
    sequence_button = st.sidebar.button("Analyze Sequence")

    if sequence_button and sequence_input:
        show_progress_bar()
        analyze_protein_sequence(sequence_input)

    # Button to trigger download preparation
    if st.sidebar.button("Fetch and Prepare Download") and global_fasta_sequence:
        st.download_button(
            label="Download FASTA",
            data=global_fasta_sequence,
            file_name=f"{protein_id}.fasta",
            mime="text/plain"
        )


# Function to fetch the protein sequence in FASTA format
@st.cache_data
def fetch_fasta(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.encode('utf-8')  # Encode as bytes for download
    else:
        st.error("Failed to retrieve the protein sequence.")
        return None


# Function to simulate a progress bar
def show_progress_bar():
    progress_text = "Operation in progress. Please wait."
    my_bar = st.progress(0, text=progress_text)

    for percent_complete in range(100):
        time.sleep(0.01)  # Simulates progress
        my_bar.progress(percent_complete + 1, text=progress_text)

    time.sleep(1)  # Pause for a moment after completion
    my_bar.empty()


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


# Display the Protein-Protein Interaction Network
def display_ppi_network(uniprot_id):
    st.subheader("Protein-Protein Interaction Network")
    G = nx.Graph()
    G.add_edge("Protein1", "Protein2")  # Placeholder for actual PPI data
    G.add_edge("Protein1", "Protein3")
    pos = nx.spring_layout(G)
    plt.figure(figsize=(8, 8))
    nx.draw(G, pos, with_labels=True, node_color='skyblue', edge_color='#FF5733', node_size=2000, font_size=10)
    st.pyplot(plt.gcf())
    plt.clf()


# Analyze a protein sequence by calculating molecular weight and alignment
def analyze_protein_sequence(sequence):
    seq = Seq(sequence)
    st.write("Molecular Weight: {:.2f} Da".format(molecular_weight(seq, seq_type='protein')))
    align_sequences(global_fasta_sequence, sequence)  # Align against the FASTA sequence

# Extract the sequence portion from the FASTA format and align two sequences
def align_sequences(fasta_seq, input_seq):
    # Debugging: Print the raw FASTA sequence for validation
    st.text("FASTA Sequence from UniProt (first 500 chars):")
    st.text(fasta_seq[:500])

    # Split the FASTA data into lines (ignoring the first header line)
    fasta_lines = fasta_seq.splitlines()
    uniprot_seq = "".join(line.strip().upper() for line in fasta_lines if not line.startswith(">"))
    input_seq = input_seq.strip().upper()

    # Debugging: Print the extracted and input sequences for validation
    st.text("Extracted UniProt Sequence (first 500 chars):")
    st.text(uniprot_seq[:500])
    st.text("Input Sequence:")
    st.text(input_seq)

    # Check for empty sequences
    if not uniprot_seq or not input_seq:
        st.warning("One or both sequences are empty. Please verify the data.")
        return

    alignments = pairwise2.align.globalxx(uniprot_seq, input_seq)
    if alignments:
        alignment_text = pairwise2.format_alignment(*alignments[0])
        st.text("Alignment:")
        st.text(alignment_text)
    else:
        st.warning("No alignments found.")





if __name__ == "__main__":
    main()
