import os
import pandas as pd
import streamlit as st
from Bio import Entrez, Medline
from io import BytesIO
import itertools
from pyvis.network import Network
import gdown
import zipfile
import spacy
import plotly.graph_objs as go

# Initialize PubMed Email
Entrez.email = st.text_input("Enter your email (for PubMed access):")

# Set up PubMed article types
article_types = {
    "Clinical Trials": "Clinical Trial[pt]",
    "Meta-Analysis": "Meta-Analysis[pt]",
    "Randomized Controlled Trials": "Randomized Controlled Trial[pt]",
    "Reviews": "Review[pt]",
    "Systematic Reviews": "Systematic Review[pt]"
}

# Helper functions
def construct_query(search_term, mesh_term, choice):
    chosen_article_type = article_types[choice]
    query = f"({search_term}) AND {chosen_article_type}"
    if mesh_term:
        query += f" AND {mesh_term}[MeSH Terms]"
    return query

def fetch_abstracts(query, num_articles):
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=num_articles)
        result = Entrez.read(handle)
        handle.close()

        ids = result['IdList']
        if not ids:
            return []

        handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
        records = Medline.parse(handle)
        articles = list(records)
        handle.close()
        return articles
    except Exception as e:
        st.write(f"Error fetching articles: {e}")
        return []

def save_to_excel(articles):
    output = BytesIO()
    data = [{
        'Title': article.get('TI', 'No title'),
        'Authors': ', '.join(article.get('AU', 'No authors')),
        'Abstract': article.get('AB', 'No abstract'),
        'Publication Date': article.get('DP', 'No date'),
        'Journal': article.get('TA', 'No journal')
    } for article in articles]

    df = pd.DataFrame(data)
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        df.to_excel(writer, index=False)
    output.seek(0)
    return output

# PubMed UI
st.title("PubMed Research Navigator & Biomedical Entity Visualizer")
search_term = st.text_input("Enter search term:")
mesh_term = st.text_input("Optional MeSH term:")
article_choice = st.selectbox("Select article type:", list(article_types.keys()))
num_articles = st.number_input("Number of articles:", min_value=1, max_value=100, value=10)

if st.button("Fetch PubMed Articles"):
    if Entrez.email and search_term:
        query = construct_query(search_term, mesh_term, article_choice)
        articles = fetch_abstracts(query, num_articles)
        
        if articles:
            excel_data = save_to_excel(articles)
            st.download_button(
                label="Download PubMed Articles",
                data=excel_data,
                file_name="pubmed_articles.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
        else:
            st.write("No articles found.")

# Biomedical Entity Visualizer
st.write("### Upload Excel for Biomedical Entity Visualization")
uploaded_file = st.file_uploader("Upload Excel file", type=["xlsx"])

@st.cache_resource
def download_and_load_model():
    base_dir = "/tmp/model/en_ner_bc5cdr_md"
    nested_model_dir = "/tmp/model/en_ner_bc5cdr_md/en_ner_bc5cdr_md-0.4.0"
    zip_path = "/tmp/model/en_ner_bc5cdr_md.zip"
    config_path = os.path.join(nested_model_dir, "config.cfg")
    download_url = "https://drive.google.com/uc?id=1kjTjVdmtLJSu7BFWMn2HMiB7eTSdmqhy"

    if not os.path.exists(nested_model_dir) or not os.path.exists(config_path):
        os.makedirs(base_dir, exist_ok=True)
        gdown.download(download_url, zip_path, quiet=False)
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall("/tmp/model")
        if not os.path.exists(config_path):
            raise FileNotFoundError("Model extraction failed.")

    return spacy.load(nested_model_dir)

# Load model and process entities
try:
    nlp = download_and_load_model()
except Exception as e:
    st.error(f"Failed to load model: {e}")
    st.stop()

def get_bc5cdr_entities(sent, entity_types):
    doc = nlp(sent)
    return [(ent.text.lower(), ent.label_) for ent in doc.ents if ent.label_ in entity_types]

def process_abstracts_from_excel(df, entity_types, allowed_relationships):
    rows = []
    entity_to_titles = {}
    for _, row in df[['Abstract', 'Title']].dropna().iterrows():
        entities = get_bc5cdr_entities(row['Abstract'], entity_types)
        for entity in entities:
            entity_to_titles.setdefault(entity[0], set()).add(row['Title'])
        for entity1, entity2 in itertools.combinations(entities, 2):
            if (entity1[1], entity2[1]) in allowed_relationships or (entity2[1], entity1[1]) in allowed_relationships:
                rows.append({'source': entity1[0], 'target': entity2[0], 'edge': f"{entity1[1]}_to_{entity2[1]}"})
    return pd.DataFrame(rows), entity_to_titles

from pyvis.network import Network

def visualize_graph_interactive(kg_df, entity_to_titles):
    """
    Visualizes the graph of entities, their relationships, and associated paper titles interactively using Pyvis,
    with specified colors for each entity type.
    """
    net = Network(height="800px", width="100%", bgcolor="#222222", font_color="white")
    net.force_atlas_2based()

    entity_colors = {'DISEASE': '#f68b24', 'CHEMICAL': '#ffffff'}
    
    # Add all nodes first
    for entity, titles in entity_to_titles.items():
        color = entity_colors.get(entity.split('_')[0], "#999999")  # Default color if not found
        title = "<br>".join(titles)
        net.add_node(entity, title=title, color=color)

    # Add edges only if nodes exist
    for _, row in kg_df.iterrows():
        source = row['source']
        target = row['target']
        
        if source in net.get_nodes() and target in net.get_nodes():
            net.add_edge(source, target, title=row['edge'])
        else:
            print(f"Warning: Edge from {source} to {target} cannot be added because one of the nodes is missing.")
    
    # Render the network as an HTML component in Streamlit
    net.show("graph.html")
    st.components.v1.html(open("graph.html", "r").read(), height=800)




