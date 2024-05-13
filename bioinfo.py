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
st.set_page_config(page_title="Protein Data Analysis", layout="wide")

# Function to set custom styles
def set_custom_styles():
    st.markdown(
        """
        <style>
            /* Add custom styles here */
            .stApp {
                background-color: #f0f2f6;
            }
            .stTextArea, .stTextInput, .stSelectbox, .stButton {
                border-radius: 10px;
                box-shadow: 0px 2px 10px rgba(0, 0, 0, 0.1);
            }
            .stButton button {
                background-color: #6d7fcc;
                color: white;
                font-weight: bold;
                border-radius: 10px;
            }
            .stButton button:hover {
                background-color: #5c6eb9;
            }
            .stTitle {
                color: #2e3f64;
                font-size: 36px;
                font-weight: bold;
                margin-bottom: 20px;
            }
            .stHeader {
                color: #2e3f64;
                font-size: 24px;
                font-weight: bold;
                margin-top: 20px;
                margin-bottom: 10px;
            }
            .stText {
                color: #2e3f64;
                font-size: 18px;
            }
            .stProgressBar div {
                background-color: #6d7fcc;
            }
            .stPlotlyChart {
                margin-top: 20px;
            }
        </style>
        """,
        unsafe_allow_html=True,
    )

# Apply custom styles
set_custom_styles()

def main():
    st.title("🔬 Protein Data Analysis App")

    st.write("## 🧬 Choose Data Source")
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
            uniProt_id = fetch_uniprot_id(sequence_input)
            if uniProt_id:
                protein_data = fetch_protein_data(uniProt_id)
                if protein_data:
                    display_protein_info(protein_data)
                    display_ppi_network(uniProt_id)

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
    st.subheader("🔍 Protein Characteristics")
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
    st.subheader("🔗 Protein-Protein Interaction Network")
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
        st.write("No interaction data available")

if __name__ == "__main__":
    main()
