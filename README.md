# Stage_three_hackbio
Reproduction of SARS-CoV-2 infection dynamics using single-cell RNA sequencing
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 1. Load Data (Assumes you have downloaded GSE166766 files)
# For the full reproduction, you would load the .mtx and .tsv files from GEO
# adata = sc.read_mtx('matrix.mtx').T
# adata.var_names = pd.read_csv('genes.tsv', header=None, sep='\t')[0]
# adata.obs = pd.read_csv('metadata.tsv', sep='\t')

# 2. Preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-') 
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# 3. Dimensionality Reduction & Clustering (Fig 1G Reproduction)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

# 4. Cell Type Identification (Markers from the paper)
marker_genes = {
    'Ciliated': ['FOXJ1', 'TPPP3'],
    'Basal': ['KRT5', 'TP63'],
    'Club': ['SCGB1A1', 'MUC5B'],
    'Goblet': ['MUC5AC', 'SPDEF'],
    'Ionocyte': ['FOXI1', 'ASCL3']
}
sc.pl.umap(adata, color=['leiden', 'FOXJ1', 'KRT5', 'SCGB1A1'], title='Cell Type Markers')

# 5. Pseudotime Analysis (Fig 3A Reproduction)
# Setting 'Basal' cells as the root because they are the progenitors
adata.uns['iroot'] = np.flatnonzero(adata.obs['leiden'] == 'Basal_Cluster_ID')[0]
sc.tl.dpt(adata)
sc.pl.umap(adata, color=['dpt_pseudotime'], title='Differentiation Pseudotime')

# 6. Viral Expression Analysis (Fig 4A/B)
# Assuming 'SARS-CoV-2' reads are in a metadata column or specific gene
sc.pl.umap(adata, color=['SARS-CoV-2_reads', 'ACE2'], cmap='Reds', title='Infection Dynamics')
