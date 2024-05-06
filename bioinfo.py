import streamlit as st
import requests
import matplotlib.pyplot as plt
import networkx as nx
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# Function to retrieve protein data from UniProt
def get_protein_data(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.txt"
    response = requests.get(url)
    protein_data = {"ID": uniprot_id}  # Initialize with the UniProt ID
    for line in response.iter_lines():
        line = line.decode("utf-8")
        if line.startswith("ID"):
            fields = line.split()
            if len(fields) >= 2:
                protein_data["Name"] = ' '.join(fields[1:])
            else:
                protein_data["Name"] = "Not available"
        elif line.startswith("SQ"):
            weight_line = next(response.iter_lines()).decode("utf-8")
            weight = weight_line.split()[-1]
            protein_data["Weight"] = weight
        elif line.startswith("DE   RecName: Full="):
            function = line.split("Full=")[1].split(";")[0]
            protein_data["Function"] = function
        elif line.startswith("DR   SUPFAM"):
            structure = line.split(";")[1]
            protein_data["Structure"] = structure
        elif line.startswith("FT   MOD_RES"):
            length = int(line.split()[2])
            protein_data["Length"] = length
        elif line.startswith("CC   -!- SUBCELLULAR LOCATION:"):
            subcellular_location = line.split("CC   -!- SUBCELLULAR LOCATION:")[1].strip()
            protein_data["Subcellular Location"] = subcellular_location
        elif line.startswith("FT   MOD_RES"):
            ptms = line.split(";")[1:]
            protein_data["PTMs"] = [ptm.strip() for ptm in ptms]
        elif line.startswith("DR   Reactome"):
            pathway = line.split(";")[1]
            protein_data["Pathway"] = pathway
        elif line.startswith("DR   MIM"):
            disease = line.split(";")[1]
            protein_data["Disease"] = disease
    return protein_data

# Function to display protein characteristics
def display_protein_characteristics(protein_data):
    st.subheader("Protein Characteristics")
    st.write(f"UniProt ID: {protein_data['ID']}")
    st.write(f"Name: {protein_data['Name']}")
    if 'Length' in protein_data:
        st.write(f"Length: {protein_data['Length']}")
    else:
        st.write("Length information not available")
    st.write(f"Molecular Weight: {protein_data['Weight']}")
    st.write(f"Function: {protein_data['Function']}")
    if 'Structure' in protein_data:
        st.write(f"Structure: {protein_data['Structure']}")
    if 'Subcellular Location' in protein_data:
        st.write(f"Subcellular Location: {protein_data['Subcellular Location']}")
    if 'PTMs' in protein_data:
        st.write("Post-Translational Modifications (PTMs):")
        for ptm in protein_data['PTMs']:
            st.write(ptm)
    if 'Pathway' in protein_data:
        st.write(f"Pathway Involvement: {protein_data['Pathway']}")
    if 'Disease' in protein_data:
        st.write(f"Disease Association: {protein_data['Disease']}")

import requests

# Function to retrieve protein-protein interaction network from STRING DB
def get_ppi_network(uniprot_id):
    proteins = ["ProteinA", "ProteinB", "ProteinC", "ProteinD", "ProteinE", "ProteinF", "ProteinG", "ProteinH", "ProteinI", "ProteinJ"]
    interactions = []
    for i in range(len(proteins)):
        for j in range(i + 1, len(proteins)):
            interactions.append((proteins[i], proteins[j]))
    return interactions

# Function to display protein-protein interaction network
def display_ppi_network(interactions):
    st.subheader("Protein-Protein Interaction Network")
    fig, ax = plt.subplots()
    G = nx.Graph()
    G.add_edges_from(interactions)
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True, node_color='skyblue', node_size=1500, font_size=10, font_color='black', edge_color='black', linewidths=1, arrows=False)
    plt.title("Protein-Protein Interaction Network")
    st.pyplot(fig)

# Function to perform sequence alignment
def perform_sequence_alignment(input_sequence, reference_sequence):
    alignments = pairwise2.align.globalxx(input_sequence, reference_sequence, one_alignment_only=True)
    return alignments[0]

# Main function
def main():
    st.title("Protein Data Explorer")
    
    # Dropdown to choose input type
    input_type = st.selectbox("Choose Input Type", ["UniProt ID", "Protein Sequence"])
    
    if input_type == "UniProt ID":
        uniprot_id = st.text_input("Enter UniProt ID:")
        if st.button("Get Protein Data"):
            if uniprot_id:
                protein_data = get_protein_data(uniprot_id)
                display_protein_characteristics(protein_data)
                ppi_network = get_ppi_network(uniprot_id)  # Dummy implementation, replace with actual retrieval
                display_ppi_network(ppi_network)
    
    elif input_type == "Protein Sequence":
        # Protein sequence input
        protein_sequence = st.text_area("Enter Protein Sequence:")
        reference_sequence = "ACGT..."  # Reference sequence for alignment

        if st.button("Perform Sequence Alignment"):
            if protein_sequence:
                alignment = perform_sequence_alignment(protein_sequence, reference_sequence)
                st.text("Sequence Alignment:")
                st.text(format_alignment(*alignment))

if __name__ == "__main__":
    main()
