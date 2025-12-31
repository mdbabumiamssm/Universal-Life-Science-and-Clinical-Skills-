import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Add skills directory to path
sys.path.append('Skills/Genomics/Single_Cell_RNA_QC')
from qc_core import calculate_qc_metrics, detect_outliers_mad, apply_hard_threshold, filter_cells, filter_genes

# Set plotting settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Data paths
data_dir = 'tests/scRNAsedata'
samples = ['GSM3901485_BM1', 'GSM3901486_BM2', 'GSM3901487_BM3']

adatas = []

print("Starting analysis...")

# 1. Load and QC each sample
for sample in samples:
    print(f"Processing {sample}...")
    sample_path = os.path.join(data_dir, sample)
    
    # Read 10x data
    try:
        adata = sc.read_10x_mtx(sample_path, var_names='gene_symbols', cache=False)
    except Exception as e:
        print(f"Error reading {sample}: {e}")
        continue
        
    adata.var_names_make_unique()
    # Note: When concatenating, we will use the keys to label batches, 
    # but setting it here is also fine for clarity
    adata.obs['batch'] = sample
    
    # Calculate QC metrics
    # Using defaults from qc_analysis.py: 
    # mad_counts=5, mad_genes=5, mad_mt=3, mt_threshold=8, min_cells=20
    calculate_qc_metrics(adata, inplace=True)
    
    # MAD-based outlier detection
    outlier_counts = detect_outliers_mad(adata, 'total_counts', n_mads=5)
    outlier_genes = detect_outliers_mad(adata, 'n_genes_by_counts', n_mads=5)
    outlier_mt = detect_outliers_mad(adata, 'pct_counts_mt', n_mads=3)
    
    # Hard threshold for MT
    high_mt = apply_hard_threshold(adata, 'pct_counts_mt', threshold=8, operator='>')
    
    # Combine filters
    pass_qc = ~(outlier_counts | outlier_genes | outlier_mt | high_mt)
    
    print(f"  Cells passing QC: {pass_qc.sum()}/{adata.n_obs} ({pass_qc.sum()/adata.n_obs*100:.1f}%)")
    
    # Filter cells
    adata = filter_cells(adata, pass_qc)
    
    # Filter genes
    filter_genes(adata, min_cells=20, inplace=True)
    
    adatas.append(adata)

if not adatas:
    print("No data loaded. Exiting.")
    sys.exit(1)

# 2. Concatenate
print("Concatenating datasets...")
# Using label='batch' creates a column 'batch' with the keys
adata_integrated = ad.concat(adatas, label='batch', keys=samples) 

# 3. Preprocessing
print("Preprocessing...")
sc.pp.normalize_total(adata_integrated, target_sum=1e4)
sc.pp.log1p(adata_integrated)
sc.pp.highly_variable_genes(adata_integrated, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_integrated.raw = adata_integrated # Save raw data
adata_integrated = adata_integrated[:, adata_integrated.var.highly_variable]
sc.pp.scale(adata_integrated, max_value=10)
sc.tl.pca(adata_integrated, svd_solver='arpack')

# 4. Integration
print("Running integration...")
use_rep = 'X_pca'
try:
    # Try importing harmony
    # Usually sc.external.pp.harmony_integrate(adata, 'batch')
    print("Attempting Harmony integration...")
    # sc.external.pp.harmony_integrate is available in newer scanpy versions
    # Sometimes needs 'harmony-pytorch' or 'harmony' python package installed.
    # We will try it, if it fails, fallback to standard PCA
    sc.external.pp.harmony_integrate(adata_integrated, 'batch')
    use_rep = 'X_pca_harmony'
    print("Harmony integration successful.")
except Exception as e:
    print(f"Harmony integration skipped (using PCA): {e}")
    use_rep = 'X_pca'

# 5. UMAP
print(f"Running UMAP using {use_rep}...")
sc.pp.neighbors(adata_integrated, n_neighbors=10, n_pcs=40, use_rep=use_rep)
sc.tl.umap(adata_integrated)

# 6. Check for CD34
print("Checking for CD34...")
genes_to_plot = ['batch']
if 'CD34' in adata_integrated.raw.var_names:
    print("CD34 found in dataset.")
    genes_to_plot.append('CD34')
else:
    print("CD34 NOT found in variable names.")

# 7. Save results
print("Saving results...")
output_dir = 'analysis_results'
os.makedirs(output_dir, exist_ok=True)

# Plotting
# sc.pl.umap saves to 'figures/' by default with 'show=False'
sc.pl.umap(adata_integrated, color=genes_to_plot, save='_integrated_CD34.png', show=False)

# Move the saved plot to output_dir
fig_path = 'figures/umap_integrated_CD34.png'
if os.path.exists(fig_path):
    os.rename(fig_path, os.path.join(output_dir, 'final_integrated_umap_CD34.png'))
    print(f"UMAP saved to {os.path.join(output_dir, 'final_integrated_umap_CD34.png')}")

adata_integrated.write(os.path.join(output_dir, 'integrated_analyzed.h5ad'))
print(f"Anndata saved to {os.path.join(output_dir, 'integrated_analyzed.h5ad')}")

print("Process complete.")
