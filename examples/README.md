# scio Examples (v0.1.2)

This directory contains example scripts demonstrating how to use the scio package for cross-platform single-cell data exchange between Python and R.

**v0.1.2**: Uses binary CSC format - 3.4x faster than H5AD in R, 14.2x faster than MTX.

## Files

| File | Description |
|------|-------------|
| `01_create_sample_data.py` | Creates synthetic single-cell data and saves as `.scio` and `.h5ad` |
| `02_python_usage.py` | Python usage examples (read, write, incremental update) |
| `03_r_usage.R` | R usage examples (Seurat, SingleCellExperiment) |

## Quick Start

### 1. Create sample data

```bash
cd examples
python 01_create_sample_data.py
```

This creates:
- `sample_data.scio/` - scio format (folder)
- `sample_data.h5ad` - H5AD format (for comparison)

### 2. Run Python examples

```bash
python 02_python_usage.py
```

### 3. Run R examples

```bash
Rscript 03_r_usage.R
```

## Cross-Platform Workflow

### Python to R

```python
# Python: Save AnnData
import scio
import scanpy as sc

adata = sc.read_h5ad("data.h5ad")
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

scio.write(adata, "processed_data.scio")
```

```r
# R: Load as Seurat
library(scio)

seurat <- scio_read("processed_data.scio", output = "Seurat")
# All embeddings (PCA, UMAP) are preserved!
DimPlot(seurat, reduction = "umap")
```

### R to Python

```r
# R: Save Seurat
library(scio)
library(Seurat)

seurat <- CreateSeuratObject(counts = pbmc.data)
seurat <- NormalizeData(seurat)
seurat <- RunPCA(seurat)
seurat <- RunUMAP(seurat, dims = 1:30)

scio_write(seurat, "seurat_data.scio")
```

```python
# Python: Load as AnnData
import scio

adata = scio.read("seurat_data.scio")
# All reductions are in adata.obsm
print(adata.obsm.keys())  # ['X_pca', 'X_umap']
```

## Incremental Updates

For iterative workflows, use `update=True` to only rewrite changed components:

```python
# Python
scio.write(adata, "data.scio")  # Initial save

# ... modify only obs ...
adata.obs['new_annotation'] = values

scio.write(adata, "data.scio", update=True)  # Only rewrites obs
```

```r
# R
scio_write(seurat, "data.scio")  # Initial save

# ... modify metadata ...
seurat$new_annotation <- values

scio_write(seurat, "data.scio", update = TRUE)  # Only rewrites changed parts
```

## Performance (v0.1.2)

Benchmark results for 100K cells × 36K genes:

| Format | R Read Time | vs scio |
|--------|-------------|---------|
| **scio v0.1.2** | **29.4s** | 1x |
| H5AD (zellkonverter) | 99.3s | 3.4x slower |
| MTX | 418.4s | 14.2x slower |

scio v0.1.2 uses binary CSC format with numpy arrays and MatrixExtra::t_shallow() for zero-copy transpose in R.
