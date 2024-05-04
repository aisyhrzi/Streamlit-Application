import streamlit as st
import requests
import json

def main():
    st.title("Protein Data Analysis App")
    
    # Input for UniProt ID
    protein_id = st.sidebar.text_input("Enter UniProt ID", value="P04637")  # Default ID for TP53 human
    analyze_button = st.sidebar.button("Analyze Protein")
    
    if analyze_button:
        protein_data = fetch_protein_data(protein_id)
        if protein_data:
            display_protein_info(protein_data)
    
    # Additional feature: Protein sequence analysis
    st.sidebar.subheader("Protein Sequence Analysis")
    sequence_input = st.sidebar.text_area("Enter Protein Sequence", value="")
    sequence_button = st.sidebar.button("Analyze Sequence")
    
    if sequence_button and sequence_input:
        analyze_protein_sequence(sequence_input)

def fetch_protein_data(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        response = requests.get(url)
        response.raise_for_status()  # Will raise an exception for HTTP codes 400 or 500
        data = response.json()
        return data
    except requests.exceptions.HTTPError as errh:
        st.error(f"HTTP Error: {errh.response.status_code} - {errh.response.text}")
    except requests.exceptions.ConnectionError as errc:
        st.error(f"Error Connecting: {errc}")
    except requests.exceptions.Timeout as errt:
        st.error(f"Timeout Error: {errt}")
    except requests.exceptions.RequestException as err:
        st.error(f"Unexpected Error: {err}")
    return None

def display_protein_info(data):
    st.subheader("Protein Characteristics")
    if "entry" in data and "protein" in data["entry"] and "recommendedName" in data["entry"]["protein"] and "fullName" in data["entry"]["protein"]["recommendedName"]:
        description = data["entry"]["protein"]["recommendedName"]["fullName"]["value"]
        st.write("Description:", description)
    if "sequence" in data["entry"]:
        sequence = data["entry"]["sequence"]["sequence"]
        st.write("Protein Length:", len(sequence))
    else:
        st.error("Protein data is incomplete or missing. Check API response format.")

def analyze_protein_sequence(sequence):
    # Placeholder for sequence analysis logic
    st.write("Sequence provided:", sequence)
    # Implement actual analysis here

if __name__ == "__main__":
    main()
