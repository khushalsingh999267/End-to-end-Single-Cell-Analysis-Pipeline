# app.py
import streamlit as st
import scanpy as sc
from PIL import Image

st.set_page_config(layout="wide")
st.title("🔬 Single-Cell Analysis Dashboard")

st.sidebar.header("Analysis Parameters")
st.sidebar.write("This dashboard visualizes the output of the Nextflow pipeline.")

# --- Main Panel ---
st.header("Cell Clustering Results")

# Load the analysis results
try:
    adata = sc.read("results/deep_analysis/downstream_analysis_results.h5ad")
    cluster_img = Image.open("results/deep_analysis/_clusters.png")

    st.image(cluster_img, caption="UMAP projection of cell clusters")

    st.info(f"""
    **Analysis Summary:**
    - **Number of cells:** {adata.n_obs}
    - **Number of genes:** {adata.n_vars}
    - **Detected clusters (Leiden):** {len(adata.obs['leiden'].unique())}
    """)

    # Display the raw data
    if st.checkbox("Show Raw Data (first 100 cells)"):
        st.dataframe(adata.to_df().head(100))

except FileNotFoundError:
    st.error("Results file not found. Please run the Nextflow pipeline first: `nextflow run main.nf`")