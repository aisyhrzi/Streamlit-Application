import streamlit as st
import requests
import time
import io
import difflib
import networkx as nx
from Bio import SeqIO
import matplotlib.pyplot as plt
from Bio.SeqUtils import molecular_weight
from xml.etree import ElementTree as ET

# Streamlit Page Config
st.set_page_config(page_title="Protein Data Analysis", layout="centered")

def main():
    st.title("Protein Data Analysis App")
    
    st.write("## Choose Data Source")
    data_source = st.radio("Select Data Source", ('UniProt ID', 'Protein Sequence'))

    if data_source == 'UniProt ID':
        protein_id = st.text_input("Enter UniProt ID", value="P04637")  # Default ID for TP53 human
        analyze_button = st.button("Analyze Protein")
        if analyze_button:
            show_progress_bar()
            protein_data = fetch_protein_data(protein_id)
            if protein_data:
                display_protein_info(protein_data)
                display_ppi_network(protein_id)
    elif data_source == 'Protein Sequence':
         sequence_input = st.text_area("Enter Protein Sequence", value="")
         analyze_button = st.button("Analyze Sequence")
         if analyze_button and sequence_input:
            show_progress_bar()
            analyze_protein_sequence(sequence_input)

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

# Function to analyze the protein sequence
def analyze_protein_sequence(sequence_input):
    # Parse the input sequence
    try:
        # Create a temporary file-like object to mimic a file handle
        with io.StringIO(sequence_input) as handle:
            records = list(SeqIO.parse(handle, "fasta"))
        
        if len(records) != 1:
            st.error("Please provide a single protein sequence in FASTA format.")
            return
        protein_sequence = str(records[0].seq)
    except Exception as e:
        st.error("Error parsing the input sequence: {}".format(str(e)))
        return

    # Calculate molecular weight
    weight = molecular_weight(protein_sequence, seq_type='protein')

    # Display results
    st.subheader("Analysis Results")
    st.write("Protein Length:", len(protein_sequence))
    st.write("Molecular Weight: {:.2f} Da".format(weight))

    # Display Protein-Protein Interaction Network
    st.subheader("Protein-Protein Interaction Network")
    ppi_data = fetch_string_ppi(protein_sequence)
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
        
if __name__ == "__main__":
    main()
