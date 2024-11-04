import pandas as pd
import spacy
import itertools
from pyvis.network import Network
import os
import gdown
import zipfile
import streamlit as st

# Function to download and load the SpaCy model from Google Drive
@st.cache_resource
def download_and_load_model():
    # Define paths
    model_dir = "/tmp/model/en_ner_bc5cdr_md"
    zip_path = "/tmp/model/en_ner_bc5cdr_md.zip"
    
    # Direct download link to the model ZIP file on Google Drive
    download_url = "https://drive.google.com/uc?id=1kjTjVdmtLJSu7BFWMn2HMiB7eTSdmqhy"

    # Download the zip file if the model folder does not exist
    if not os.path.exists(model_dir):
        os.makedirs("/tmp/model", exist_ok=True)
        gdown.download(download_url, zip_path, quiet=False)
        
        # Unzip the downloaded model folder
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall("/tmp/model")
    
    # Load the model
    return spacy.load(model_dir)

# Initialize the SpaCy model
nlp = download_and_load_model()

def get_bc5cdr_entities(sent, entity_types):
    """
    Extracts specified types of entities using the BC5CDR model and normalizes them to lowercase.
    """
    doc = nlp(sent)
    entities = [(ent.text.lower(), ent.label_) for ent in doc.ents if ent.label_ in entity_types]
    return entities

def process_abstracts_from_excel(excel_path, entity_types, allowed_relationships):
    """
    Processes each abstract to identify specified entities, their relationships, and associated paper titles,
    filtering by allowed relationships.
    """
    df = pd.read_excel(excel_path)
    rows = []
    entity_to_titles = {}  # Dictionary to map entities to their associated paper titles
    
    for _, row in df[['Abstract', 'Title']].dropna().iterrows():
        abstract, title = row['Abstract'], row['Title']
        entities = get_bc5cdr_entities(abstract, entity_types)
        for entity in entities:
            entity_to_titles.setdefault(entity[0], set()).add(title)
        for entity1, entity2 in itertools.combinations(entities, 2):
            if (entity1[1], entity2[1]) in allowed_relationships or (entity2[1], entity1[1]) in allowed_relationships:
                rows.append({'source': entity1[0], 'target': entity2[0], 'edge': f"{entity1[1]}_to_{entity2[1]}"})
    
    return pd.DataFrame(rows), entity_to_titles

def visualize_graph_interactive(kg_df, entity_to_titles, output_file):
    """
    Visualizes the graph of entities, their relationships, and associated paper titles interactively using Pyvis,
    with specified colors for each entity type.
    """
    net = Network(notebook=False, height="100%", width="100%", bgcolor="#222222", font_color="white")
    net.force_atlas_2based()
    
    entity_colors = {'DISEASE': '#f68b24', 'GENE': '#4caf50', 'CHEMICAL': '#ffffff', 'PROTEIN': '#f44336'}
    
    for _, row in kg_df.iterrows():
        source_title = "<br>".join(entity_to_titles[row['source']])
        target_title = "<br>".join(entity_to_titles[row['target']])
        # Determine the color based on entity type
        source_color = entity_colors.get(row['edge'].split('_to_')[0], "#999999")  # Default color if not found
        target_color = entity_colors.get(row['edge'].split('_to_')[1], "#999999")  # Default color if not found
        net.add_node(row['source'], title=source_title, color=source_color)
        net.add_node(row['target'], title=target_title, color=target_color)
        net.add_edge(row['source'], row['target'], title=row['edge'])
    
    net.write_html(output_file)

if __name__ == "__main__":
    excel_path = input("Enter the path to the Excel file containing abstracts: ")
    print("Enter the entity types for which you want to extract relationships (e.g., CHEMICAL, DISEASE):")
    entity_types = input("Entity types (separated by comma): ").strip().split(", ")
    
    print("Enter allowed relationships (e.g., CHEMICAL-DISEASE, DISEASE-GENE):")
    allowed_rel_input = input("Allowed relationships (separated by comma): ").strip().split(", ")
    allowed_relationships = [tuple(pair.split("-")) for pair in allowed_rel_input]
    
    kg_df, entity_to_titles = process_abstracts_from_excel(excel_path, entity_types, allowed_relationships)
    
    print(f"Processed {len(kg_df)} relationships for visualization.")
    
    output_file = input("Enter the full path and file name for the output HTML file (e.g., /path/to/output.html): ")
    visualize_graph_interactive(kg_df, entity_to_titles, output_file)


