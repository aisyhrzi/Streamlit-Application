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
st.set_page_config(page_title="Enhanced Protein Data Analysis", layout="wide")

# Main application function
def main():
    st.title("Enhanced Protein Data Analysis App")
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
        analyze_protein_sequence(protein_id, sequence_input)

    # Button to trigger download preparation
    if st.sidebar.button("Fetch and Prepare Download"):
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
@st.cache
def fetch_fasta(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
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

# Function to fetch protein data and parse XML with detailed fields
def fetch_protein_data(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
    response = requests.get(url)
    if response.status_code == 200:
        root = ET.fromstring(response.content)
        description = root.findtext('.//{http://uniprot.org/uniprot}fullName')
        sequence = root.findtext('.//{http://uniprot.org/uniprot}sequence').replace('\n', '').strip()
        organism = root.findtext('.//{http://uniprot.org/uniprot}organism/{http://uniprot.org/uniprot}name[@type="scientific"]')
        function = root.findtext('.//{http://uniprot.org/uniprot}comment[@type="function"]/{http://uniprot.org/uniprot}text')
        subcellular_location = root.findtext('.//{http://uniprot.org/uniprot}comment[@type="subcellular location"]/{http://uniprot.org/uniprot}location')
        pathway = root.findtext('.//{http://uniprot.org/uniprot}dbReference[@type="Reactome"]')
        disease = root.findtext('.//{http://uniprot.org/uniprot}comment[@type="disease"]/{http://uniprot.org/uniprot}text')

        return {
            "description": description,
            "sequence": sequence,
            "organism": organism,
            "function": function,
            "subcellular_location": subcellular_location,
            "pathway": pathway,
            "disease": disease
        }
    else:
        st.error("Failed to retrieve data.")
        return None
    # Display protein characteristics
def display_protein_info(data):
    st.subheader("Protein Characteristics")
    st.write("Description:", data["description"])
    st.write("Organism:", data["organism"])
    st.write("Protein Length:", len(data["sequence"]))
    st.write("Molecular Weight: {:.2f} Da".format(molecular_weight(data["sequence"], seq_type='protein')))
    st.write("Function:", data["function"])
    st.write("Subcellular Location:", data["subcellular_location"])
    st.write("Pathway:", data["pathway"])
    st.write("Disease Association:", data["disease"])

# Display the Protein-Protein Interaction Network
def display_ppi_network(uniprot_id):
    st.subheader("Protein-Protein Interaction Network")
    ppi_data = fetch_string_ppi(uniprot_id)
    if ppi_data:
        G = nx.Graph()
        for interaction in ppi_data:
            protein1 = interaction["preferredName_A"]
            protein2 = interaction["preferredName_B"]
            score = interaction["score"]
            G.add_edge(protein1, protein2, weight=score)

        pos = nx.spring_layout(G, k=0.5)  # Adjust k to change spacing between nodes
        plt.figure(figsize=(10, 10))
        nx.draw(G, pos, with_labels=True, node_color='skyblue', edge_color='#FF5733', node_size=2000, font_size=10, width=[data['weight'] for _, _, data in G.edges(data=True)])
        st.pyplot(plt.gcf())
        plt.clf()
    else:
        st.write("No interaction data available.")

# Function to fetch protein-protein interaction data from STRING DB
def fetch_string_ppi(uniprot_id, min_score=700):
    url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": uniprot_id,  # Protein identifier
        "species": 9606,            # Species (9606 for Homo sapiens)
        "required_score": min_score  # Minimum interaction score
    }
    response = requests.get(url, params=params)
    if response.status_code == 200:
        return response.json()
    else:
        st.error("Failed to retrieve data from STRING.")
        return None

# Analyze a protein sequence by calculating molecular weight and alignment
def analyze_protein_sequence(protein_id, sequence):
    fasta_sequence = fetch_fasta(protein_id)
    if fasta_sequence:
        seq1 = Seq(fasta_sequence.decode("utf-8").split('\n', 1)[1].replace('\n', ''))
        seq2 = Seq(sequence)
        st.write("Molecular Weight: {:.2f} Da".format(molecular_weight(seq2, seq_type='protein')))
        align_sequences(seq1, seq2) 
    else:
        st.error("Failed to fetch protein sequence.")

# Align two sequences and show the alignment
def align_sequences(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    if alignments:
        alignment_text = pairwise2.format_alignment(*alignments[0])
        st.text("Alignment:")
        st.text(alignment_text)
    else:
        st.write("No alignment found.")

if name == "main":
    main()
