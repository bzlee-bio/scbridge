#!/usr/bin/env python3
"""
Python usage examples for scio

Demonstrates how to use scio to save and load AnnData objects.
"""
import scio

# =============================================================================
# Basic Usage
# =============================================================================

print("=" * 60)
print("scio Python Usage Examples")
print("=" * 60)

# Load data from .scio format
print("\n1. Loading data from .scio format...")
adata = scio.read("sample_data.scio")
print(f"   Loaded: {adata.n_obs} cells x {adata.n_vars} genes")
print(f"   Embeddings: {list(adata.obsm.keys())}")
print(f"   Layers: {list(adata.layers.keys())}")

# =============================================================================
# Save with different options
# =============================================================================

print("\n2. Saving data...")

# Basic save
scio.write(adata, "output_basic.scio", overwrite=True)
print("   Saved: output_basic.scio")

# =============================================================================
# Incremental Updates
# =============================================================================

print("\n3. Incremental update example...")

# First save
scio.write(adata, "output_incremental.scio", overwrite=True)
print("   Initial save complete")

# Modify only cell metadata
adata.obs['new_column'] = 'test_value'

# Update - only changed components will be rewritten
scio.write(adata, "output_incremental.scio", update=True)
print("   Incremental update complete (only obs was rewritten)")

# =============================================================================
# Cross-platform compatibility
# =============================================================================

print("\n4. Cross-platform compatibility...")
print("   Files saved by Python can be read in R:")
print("   R> library(scio)")
print("   R> seurat <- scio_read('output_basic.scio', output='Seurat')")
print("   R> sce <- scio_read('output_basic.scio', output='SCE')")

# =============================================================================
# Load from H5AD and convert
# =============================================================================

print("\n5. Converting from H5AD to scio...")
import anndata as ad
adata_h5ad = ad.read_h5ad("sample_data.h5ad")
scio.write(adata_h5ad, "converted_from_h5ad.scio", overwrite=True)
print("   Converted: sample_data.h5ad -> converted_from_h5ad.scio")

# =============================================================================
# Cleanup
# =============================================================================
import shutil
for path in ["output_basic.scio", "output_incremental.scio", "converted_from_h5ad.scio"]:
    try:
        shutil.rmtree(path)
    except:
        pass

print("\n" + "=" * 60)
print("Examples completed successfully!")
print("=" * 60)
