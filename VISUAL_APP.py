import os
import pandas as pd
import spacy
import itertools
from pyvis.network import Network
from Bio import Entrez, Medline
from io import BytesIO
import streamlit as st
import gdown
import zipfile
import subprocess




st.write("Installed packages:")
st.write(subprocess.run(["pip", "freeze"], capture_output=True, text=True).stdout)

# Debugging Info
st.write("Python version:", os.sys.version)
st.write("Installed packages:")
st.write(subprocess.run(["pip", "freeze"], capture_output=True, text=True).stdout)




def download_and_load_model():
    import gdown
    import zipfile

    base_dir = "/tmp/model/en_ner_bc5cdr_md"
    nested_model_dir = "/tmp/model/en_ner_bc5cdr_md/en_ner_bc5cdr_md-0.4.0"
    zip_path = "/tmp/model/en_ner_bc5cdr_md.zip"
    config_path = os.path.join(nested_model_dir, "config.cfg")
    download_url = "https://drive.google.com/uc?id=1kjTjVdmtLJSu7BFWMn2HMiB7eTSdmqhy"

    try:
        if not os.path.exists(config_path):
            st.write("Model not found. Starting download...")
            os.makedirs(base_dir, exist_ok=True)

            output_path = gdown.download(download_url, zip_path, quiet=False)
            if not output_path or not os.path.exists(zip_path):
                raise RuntimeError("Download failed or file not found after download.")

            st.write("Download complete. Extracting...")
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall("/tmp/model")

            if not os.path.exists(config_path):
                raise FileNotFoundError("Extraction complete, but config.cfg not found.")

        st.write("Loading SpaCy model...")
        return spacy.load(nested_model_dir)

    except Exception as e:
        st.error(f"Failed to download or load the model: {e}")
        raise e


# PubMed Search Functions
def construct_query(search_term, mesh_term, choice):
    article_types = {
        "Clinical Trial": "Clinical Trial[pt]",
        "Meta-Analysis": "Meta-Analysis[pt]",
        "Randomized Controlled Trial": "Randomized Controlled Trial[pt]",
        "Review": "Review[pt]",
        "Case Reports": "Case Reports[pt]",
        "Comparative Study": "Comparative Study[pt]",
        "Observational Study": "Observational Study[pt]",
        "Controlled Clinical Trial": "Controlled Clinical Trial[pt]",
        "Evaluation Study": "Evaluation Study[pt]",
        "Multicenter Study": "Multicenter Study[pt]",
        "Practice Guideline": "Practice Guideline[pt]",
        "Validation Study": "Validation Study[pt]"
    }

    query = f"({search_term})"

    if choice != "All Articles" and choice in article_types:
        query += f" AND {article_types[choice]}"

    if mesh_term:
        query += f" AND {mesh_term}[MeSH Terms]"

    return query



def fetch_abstracts(query, num_articles, email):
    Entrez.email = email
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
    return output, df

# Entity Extraction and Visualization Functions
def get_bc5cdr_entities(sent, entity_types):
    doc = nlp(sent)
    return [(ent.text.lower(), ent.label_) for ent in doc.ents if ent.label_ in entity_types]

def process_abstracts_from_excel(df, entity_types, allowed_relationships):
    rows = []
    entity_to_titles = {}
    for _, row in df[['Abstract', 'Title']].dropna().iterrows():
        entities = get_bc5cdr_entities(row['Abstract'], entity_types)
        for entity_text, entity_type in entities:
            # Store each entity with its type and associated titles
            if entity_text not in entity_to_titles:
                entity_to_titles[entity_text] = {"titles": set(), "type": entity_type}
            entity_to_titles[entity_text]["titles"].add(row['Title'])
        
        # Generate relationships based on allowed types
        for entity1, entity2 in itertools.combinations(entities, 2):
            if (entity1[1], entity2[1]) in allowed_relationships or (entity2[1], entity1[1]) in allowed_relationships:
                rows.append({
                    'source': entity1[0], 'target': entity2[0], 
                    'edge': f"{entity1[1]}_to_{entity2[1]}"
                })
    
    return pd.DataFrame(rows), entity_to_titles

def visualize_graph_interactive(kg_df, entity_to_titles):
    from pyvis.network import Network

    net = Network(height="100vh", width="100vw", bgcolor="#222222", font_color="white")

    entity_colors = {
        "CHEMICAL": "green",
        "DISEASE": "red"
    }

    for entity, details in entity_to_titles.items():
        entity_type = details["type"]
        color = entity_colors.get(entity_type.upper(), "#999999")
        title = "<br>".join(details["titles"])
        net.add_node(entity, title=title, color=color)

    for _, row in kg_df.iterrows():
        if row['source'] in net.get_nodes() and row['target'] in net.get_nodes():
            net.add_edge(row['source'], row['target'], title=row['edge'])

    net.force_atlas_2based(
        gravity=-60,
        central_gravity=0.002,
        spring_length=100,
        spring_strength=0.01,
        damping=0.6
    )

    html_path = "graph_download.html"
    net.save_graph(html_path)

    # Read the saved HTML
    with open(html_path, "r") as file:
        html_content = file.read()

    # Inject JavaScript for click-based filtering
    


    html_content = html_content.replace("</body>", custom_js + "\n</body>")

    # Add full-page CSS styling
    full_page_css = """
    <style>
        body, html {
            margin: 0;
            padding: 0;
            width: 100vw;
            height: 100vh;
            overflow: hidden;
        }
        #mynetwork {
            width: 100vw;
            height: 100vh;
        }
    </style>
    """
    html_content = html_content.replace("<head>", f"<head>{full_page_css}")

    # Save modified HTML
    with open(html_path, "w") as file:
        file.write(html_content)

    return html_path


# Streamlit UI
st.title("PubMed Research Navigator & Biomedical Entity Visualizer")
# Ensure SpaCy model is loaded and available globally via session state
if "nlp_model" not in st.session_state:
    try:
        st.session_state["nlp_model"] = download_and_load_model()
        st.write("BIOMEDICAL model loaded and stored.")
    except Exception as e:
        st.error(f"Failed to load model: {e}")
        st.stop()

# Assign globally used nlp
nlp = st.session_state["nlp_model"]

# Store data in session state to persist across interactions
if "df" not in st.session_state:
    st.session_state["df"] = None
if "html_path" not in st.session_state:
    st.session_state["html_path"] = None

# PubMed Input Fields
email = st.text_input("Enter your email for PubMed access:")
search_term = st.text_input("Enter search term:")
mesh_term = st.text_input("Optional MeSH term:")
article_choice = st.selectbox(
    "Select article type:",
    [
        "All Articles",
        "Clinical Trial",
        "Meta-Analysis",
        "Randomized Controlled Trial",
        "Review",
        "Case Reports",
        "Comparative Study",
        "Observational Study",
        "Controlled Clinical Trial",
        "Evaluation Study",
        "Multicenter Study",
        "Practice Guideline",
        "Validation Study"
    ]
)

num_articles = st.number_input("Number of articles to fetch:", min_value=1, max_value=100, value=10)

# Fetch Articles and Download Excel
if st.button("Fetch PubMed Articles"):
    if email and search_term:
        query = construct_query(search_term, mesh_term, article_choice)
        articles = fetch_abstracts(query, num_articles, email)
        if articles:
            excel_data, df = save_to_excel(articles)
            st.session_state["df"] = df  # Save the DataFrame in session state
            st.session_state["excel_data"] = excel_data  # Save excel data in session state
            st.write("Excel file ready for entity extraction.")

# Show Download Button for Excel File if Data Exists
if "excel_data" in st.session_state and st.session_state["excel_data"]:
    st.download_button(
        label="Download PubMed Articles as Excel",
        data=st.session_state["excel_data"],
        file_name="pubmed_articles.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    )

# Entity Extraction and Visualization
if st.session_state["df"] is not None:
    df = st.session_state["df"]
    entity_types_input = st.text_input("Enter entity types (e.g., CHEMICAL, DISEASE)", "CHEMICAL, DISEASE")
    allowed_rel_input = st.text_input("Enter allowed relationships (e.g., CHEMICAL-DISEASE)")
    
    if st.button("Process and Generate Graph for Download"):
        entity_types = [et.strip() for et in entity_types_input.split(",")]
        allowed_relationships = [tuple(rel.strip().split("-")) for rel in allowed_rel_input.split(",") if "-" in rel]
        
        kg_df, entity_to_titles = process_abstracts_from_excel(df, entity_types, allowed_relationships)
        st.write(f"Processed {len(kg_df)} relationships for visualization.")
        
        # Generate HTML graph and save to session state
        html_path = visualize_graph_interactive(kg_df, entity_to_titles)
        st.session_state["html_path"] = html_path

# Provide download button for HTML if generated
if st.session_state["html_path"]:
    with open(st.session_state["html_path"], "r") as file:
        html_content = file.read()
        st.download_button(
            label="Download HTML Visualization",
            data=html_content,
            file_name="entity_relationship_graph.html",
            mime="text/html"
        )
    
    # Optional cleanup after download
    if os.path.exists(st.session_state["html_path"]):
        os.remove(st.session_state["html_path"])
        st.session_state["html_path"] = None


