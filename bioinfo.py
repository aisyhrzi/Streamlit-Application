import streamlit as st
import requests
import networkx as nx
import matplotlib.pyplot as plt
from Bio.SeqUtils import molecular_weight
from Bio import pairwise2
from Bio.Seq import Seq
from xml.etree import ElementTree as ET

# Streamlit Page Config
st.set_page_config(page_title="Protein Data Analysis", layout="wide")

def main():
    st.title("Protein Data Analysis App")
    protein_id = st.sidebar.text_input("Enter UniProt ID", value="P04637")  # Default ID for TP53 human

    analyze_button = st.sidebar.button("Analyze Protein")

    if analyze_button:
        protein_data = fetch_protein_data(protein_id)
        if protein_data:
            display_protein_info(protein_data)
            display_ppi_network(protein_id)

    st.sidebar.subheader("Protein Sequence Analysis")
    sequence_input = st.sidebar.text_area("Enter Protein Sequence", value="")
    sequence_button = st.sidebar.button("Analyze Sequence")

    if sequence_button and sequence_input:
        analyze_protein_sequence(sequence_input)

def fetch_protein_data(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.xml"
    response = requests.get(url)
    if response.status_code == 200:
        root = ET.fromstring(response.content)
        # Parsing XML data to retrieve specific information
        description = root.findtext('.//{http://uniprot.org/uniprot}fullName')
        sequence = root.findtext('.//{http://uniprot.org/uniprot}sequence').replace('\n', '').strip()
        return {
            "description": description,
            "sequence": sequence
        }
    else:
        st.error("Failed to retrieve data.")
        return None

def display_protein_info(data):
    st.subheader("Protein Characteristics")
    st.write("Description:", data["description"])
    st.write("Protein Length:", len(data["sequence"]))
    st.write("Molecular Weight: {:.2f} Da".format(molecular_weight(data["sequence"], seq_type='protein')))

def display_ppi_network(uniprot_id):
    st.subheader("Protein-Protein Interaction Network")
    # Replace placeholder edges with real PPI data if available
    G = nx.Graph()
    G.add_edge("Protein1", "Protein2")  # Placeholder for actual PPI data
    G.add_edge("Protein1", "Protein3")
    pos = nx.spring_layout(G)
    plt.figure(figsize=(8, 8))
    nx.draw(G, pos, with_labels=True, node_color='skyblue', edge_color='#FF5733', node_size=2000, font_size=10)
    st.pyplot(plt.gcf())
    plt.clf()

def analyze_protein_sequence(sequence):
    seq = Seq(sequence)
    st.write("Molecular Weight: {:.2f} Da".format(molecular_weight(seq, seq_type='protein')))
    align_sequences("MVMEESQTSDQSKE", sequence)  # Example alignment with dummy data

def align_sequences(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    alignment_text = pairwise2.format_alignment(*alignments[0])
    st.text("Alignment:")
    st.text(alignment_text)

if __name__ == "__main__":
    main()
