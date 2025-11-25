#!/usr/bin/env python3
"""
Create sample single-cell data for scio examples

This script creates a small synthetic dataset that demonstrates
all features of the scio format.
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scio

# Set random seed for reproducibility
np.random.seed(42)

# =============================================================================
# Create synthetic single-cell data
# =============================================================================
n_cells = 500
n_genes = 200

print("Creating sample single-cell dataset...")
print(f"  Cells: {n_cells}")
print(f"  Genes: {n_genes}")

# 1. Expression matrix (sparse)
# Simulate count data with some structure
counts = np.random.negative_binomial(n=1, p=0.1, size=(n_cells, n_genes))
counts = sp.csr_matrix(counts)

# 2. Cell metadata (obs)
cell_types = np.random.choice(['T_cell', 'B_cell', 'Monocyte', 'NK_cell'], n_cells)
samples = np.random.choice(['Sample_A', 'Sample_B', 'Sample_C'], n_cells)
n_counts = np.array(counts.sum(axis=1)).flatten()
n_genes_detected = np.array((counts > 0).sum(axis=1)).flatten()

obs = pd.DataFrame({
    'cell_type': pd.Categorical(cell_types),
    'sample': pd.Categorical(samples),
    'n_counts': n_counts,
    'n_genes': n_genes_detected,
    'percent_mito': np.random.uniform(0, 10, n_cells),
}, index=[f'Cell_{i:04d}' for i in range(n_cells)])

# 3. Gene metadata (var)
gene_names = [f'Gene_{i:03d}' for i in range(n_genes)]
var = pd.DataFrame({
    'gene_symbol': gene_names,
    'mean_expression': np.array(counts.mean(axis=0)).flatten(),
    'highly_variable': np.random.choice([True, False], n_genes, p=[0.1, 0.9]),
}, index=gene_names)

# 4. Create AnnData object
adata = ad.AnnData(
    X=counts,
    obs=obs,
    var=var,
)

# 5. Add embeddings (obsm)
# Simulate PCA
adata.obsm['X_pca'] = np.random.randn(n_cells, 50).astype(np.float32)
# Simulate UMAP
adata.obsm['X_umap'] = np.random.randn(n_cells, 2).astype(np.float32)

# 6. Add cell-cell graphs (obsp)
# Simulate nearest neighbor graph
from scipy.sparse import random as sparse_random
connectivities = sparse_random(n_cells, n_cells, density=0.01, format='csr')
connectivities = (connectivities + connectivities.T) / 2  # Make symmetric
adata.obsp['connectivities'] = connectivities

# 7. Add layers
adata.layers['normalized'] = sp.csr_matrix(np.log1p(counts.toarray()))

# 8. Add raw data
adata.raw = adata.copy()

# 9. Add unstructured metadata (uns)
adata.uns['dataset_name'] = 'sample_pbmc'
adata.uns['scio_version'] = scio.__version__
adata.uns['cell_type_colors'] = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
adata.uns['pca'] = {
    'variance_ratio': np.random.uniform(0.01, 0.1, 50).tolist(),
}

print("\nAnnData object created:")
print(adata)

# =============================================================================
# Save to different formats
# =============================================================================

# Save as .scio (our format)
print("\n" + "="*60)
print("Saving to .scio format...")
scio.write(adata, "sample_data.scio", overwrite=True)
print("Saved: sample_data.scio")

# Save as .h5ad (for comparison)
print("\nSaving to .h5ad format...")
adata.write_h5ad("sample_data.h5ad")
print("Saved: sample_data.h5ad")

print("\n" + "="*60)
print("Sample data files created successfully!")
print("\nYou can now use these files to test:")
print("  - Python: scio.read('sample_data.scio')")
print("  - R: scio_read('sample_data.scio', output='Seurat')")
print("  - R: scio_read('sample_data.scio', output='SCE')")
