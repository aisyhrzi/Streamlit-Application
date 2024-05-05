import streamlit as st
import requests
import networkx as nx
import matplotlib.pyplot as plt
from Bio.SeqUtils import molecular_weight
from Bio import pairwise2
from Bio.Seq import Seq
from xml.etree import ElementTree as ET
import time

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
            display_protein_info(protein_data)
            display_ppi_network(protein_id)

    st.sidebar.subheader("Protein Sequence Analysis")
    sequence_input = st.sidebar.text_area("Enter Protein Sequence", value="")
    sequence_button = st.sidebar.button("Analyze Sequence")

    if sequence_button and sequence_input:
        show_progress_bar()
        analyze_protein_sequence(sequence_input)

    # Button to trigger download preparation
    if st.sidebar.button("Fetch and Prepare Download"):
        show_progress_bar()  # Show progress before fetching
        fasta_data = fetch_fasta(protein_id)

        if fasta_data:
            # Download button for FASTA
            st.download_button(
                label="Download FASTA",
                data=fasta_data,
                file_name=f"{protein_id}.fasta",
                mime="text/plain"
            )
        else:
            st.error("Failed to retrieve the protein sequence.")

# Function to fetch the protein sequence in FASTA format
@st.cache_data
def fetch_fasta(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.encode('utf-8')  # Encode as bytes for download
    else:
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
    align_sequences("MVMEESQTSDQSKE", sequence)  # Example alignment with dummy data

# Align two sequences and show the alignment
def align_sequences(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    alignment_text = pairwise2.format_alignment(*alignments[0])
    st.text("Alignment:")
    st.text(alignment_text)

if __name__ == "__main__":
    main()
