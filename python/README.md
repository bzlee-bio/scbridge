# scio - Python Package

Cross-platform single-cell RNA-seq data storage for Python (AnnData) and R (Seurat/SingleCellExperiment).

## Installation

### From GitHub

```bash
pip install git+https://github.com/bzlee-bio/scio.git#subdirectory=python
```

### Local Development

```bash
cd scio/python
pip install -e .
```

## Quick Start

```python
import anndata as ad
import scio

# Load your data
adata = ad.read_h5ad("data.h5ad")

# Save to .scio format (works with R!)
scio.write(adata, "data.scio")

# Load back
adata = scio.read("data.scio")
```

## Features

- **Complete data preservation**: Saves ALL 10 AnnData components
  - X (expression matrix)
  - obs (cell metadata)
  - var (gene metadata)
  - obsm (cell embeddings: PCA, UMAP, etc.)
  - varm (gene embeddings)
  - obsp (cell-cell graphs)
  - varp (gene-gene graphs)
  - layers (additional matrices)
  - raw (raw counts)
  - uns (unstructured metadata)

- **Incremental updates**: Hash-based change detection
  - Only rewrites modified components
  - Dramatically faster for iterative workflows
  - Use `update=True` parameter

- **Cross-platform compatibility**: Works seamlessly with R
  - Load in R as Seurat or SingleCellExperiment objects
  - See `../R/README.md` for R usage

- **Efficient storage**:
  - MTX format for sparse matrices (universal compatibility)
  - Parquet format for metadata (preserves dtypes, 10-50x faster than CSV)
  - Direct folder access (no extraction needed)

- **Large dataset support**: Optimized for 1M+ cells

## API Reference

### write()

```python
scio.write(adata, path, overwrite=False, update=False, compress=True)
```

Save AnnData object to .scio folder.

**Parameters:**
- `adata` (AnnData): AnnData object to save
- `path` (str or Path): Output .scio folder path
- `overwrite` (bool): Whether to overwrite existing folder (default: False)
- `update` (bool): Whether to perform incremental update using hash-based change detection (default: False). Only changed components will be rewritten.
- `compress` (bool): Whether to gzip compress MTX matrix files (default: True). Compression reduces file size by ~3-5x but is slower to write. TSV files are always gzipped.

**Example:**
```python
scio.write(adata, "data.scio")
scio.write(adata, "data.scio", overwrite=True)

# Incremental update (only writes changed components)
scio.write(adata, "data.scio", update=True)

# Disable MTX compression for faster writes (larger files)
scio.write(adata, "data.scio", compress=False)
```

### read()

```python
adata = scio.read(path)
```

Load AnnData object from .scio folder (or legacy tar archive).

**Parameters:**
- `path` (str or Path): Path to .scio folder (also supports legacy tar archives)

**Returns:**
- `AnnData`: Reconstructed AnnData object

**Example:**
```python
adata = scio.read("data.scio")
```

## File Format

The .scio folder structure:

```
data.scio/
├── manifest.json           # Metadata about saved components
├── matrix.mtx[.gz]         # X (expression matrix, optionally gzipped)
├── barcodes.tsv.gz         # Cell IDs (always gzipped)
├── features.tsv.gz         # Gene IDs (always gzipped)
├── obs.parquet             # Cell metadata
├── var.parquet             # Gene metadata
├── obsm/                   # Cell embeddings
│   ├── X_pca.parquet
│   └── X_umap.parquet
├── varm/                   # Gene embeddings (if any)
├── obsp/                   # Cell-cell graphs
│   ├── distances.mtx[.gz]
│   └── connectivities.mtx[.gz]
├── varp/                   # Gene-gene graphs (if any)
├── layers/                 # Additional matrices (if any)
├── raw/                    # Raw data (if present)
│   ├── matrix.mtx[.gz]
│   ├── features.tsv.gz
│   └── var.parquet
└── uns.json                # Unstructured metadata
```

### Compression

The `compress` parameter controls **MTX file compression only**:

| File Type | compress=True | compress=False |
|-----------|---------------|----------------|
| MTX (matrix, graphs, layers) | `.mtx.gz` (smaller, slower) | `.mtx` (larger, faster) |
| TSV (barcodes, features) | `.tsv.gz` (always) | `.tsv.gz` (always) |
| Parquet (metadata, embeddings) | Built-in compression | Built-in compression |
| JSON (uns) | No compression | No compression |

**When to use `compress=False`:**
- Large datasets where write speed is critical
- Temporary files that will be deleted soon
- When disk space is not a concern

## Cross-Platform Usage

### Python to R

```python
# Python: Save
import anndata as ad
import scio

adata = ad.read_h5ad("data.h5ad")
scio.write(adata, "data.scio")
```

```r
# R: Load as Seurat
library(scio)

seurat_obj <- scio_read("data.scio", output = "Seurat")
```

### R to Python

```r
# R: Save from Seurat
library(scio)
library(Seurat)

scio_write(seurat_obj, "data.scio")
```

```python
# Python: Load
import scio

adata = scio.read("data.scio")
```

## Requirements

- Python >= 3.8
- anndata >= 0.8.0
- pandas >= 1.5.0
- numpy >= 1.21.0
- scipy >= 1.7.0
- pyarrow >= 10.0.0

## License

MIT License

## Contributing

See [CONTRIBUTING.md](../../CONTRIBUTING.md)

<!-- ## Citation

If you use scio in your research, please cite:

```
[Citation information will be added]
``` -->
