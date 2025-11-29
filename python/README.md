# scio - Python Package (v0.1.3)

Cross-platform single-cell RNA-seq data storage for Python (AnnData) and R (Seurat/SingleCellExperiment).

**v0.1.3 Highlights:**
- Binary CSR/CSC format for fast sparse matrix I/O
- Configurable sparse format: CSR (fastest write) or CSC (faster R read)
- 3.4x faster than H5AD in R, 14.2x faster than MTX
- cells×genes orientation eliminates transpose on Python read

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

# Or use CSC format for faster R reading
scio.write(adata, "data.scio", sparse_format='csc')

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

- **Fast binary CSR/CSC format** (v0.1.3):
  - Sparse matrices stored as numpy arrays (.npy files)
  - CSR (default): Fastest write from Python (scanpy outputs CSR)
  - CSC: Faster R read performance
  - cells×genes orientation - no transpose on Python read
  - R read: 29.4s vs H5AD 99.3s vs MTX 418.4s

- **Incremental updates**: Hash-based change detection
  - Only rewrites modified components
  - Dramatically faster for iterative workflows
  - Use `update=True` parameter

- **Cross-platform compatibility**: Works seamlessly with R
  - Load in R as Seurat or SingleCellExperiment objects
  - See `../R/README.md` for R usage

- **Efficient storage**:
  - Binary numpy arrays for sparse matrices
  - Parquet format for metadata (preserves dtypes)
  - Direct folder access (no extraction needed)

- **Large dataset support**: Optimized for 1M+ cells

## API Reference

### write()

```python
scio.write(adata, path, sparse_format='csr', overwrite=False, update=False)
```

Save AnnData object to .scio folder using binary sparse format.

**Parameters:**
- `adata` (AnnData): AnnData object to save
- `path` (str or Path): Output .scio folder path
- `sparse_format` (str): Format for sparse matrices - `'csr'` (default, fastest write) or `'csc'` (faster R read)
- `overwrite` (bool): Whether to overwrite existing folder (default: False)
- `update` (bool): Whether to perform incremental update using hash-based change detection (default: False). Only changed components will be rewritten.

**Example:**
```python
scio.write(adata, "data.scio")
scio.write(adata, "data.scio", overwrite=True)

# Use CSC format for faster R reading
scio.write(adata, "data.scio", sparse_format='csc')

# Incremental update (only writes changed components)
scio.write(adata, "data.scio", update=True)
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

scio v0.1.3 uses binary CSR/CSC sparse format:

```
data.scio/
├── manifest.json           # Metadata (format version, orientation)
├── matrix.data.npy         # Sparse data array (binary numpy)
├── matrix.indices.npy      # Sparse indices array (binary numpy)
├── matrix.indptr.npy       # Sparse indptr array (binary numpy)
├── shape.json              # Matrix dimensions and orientation
├── barcodes.tsv.gz         # Cell IDs
├── features.tsv.gz         # Gene IDs
├── obs.parquet             # Cell metadata
├── var.parquet             # Gene metadata
├── obsm/                   # Cell embeddings (parquet)
│   ├── X_pca.parquet
│   └── X_umap.parquet
├── varm/                   # Gene embeddings (if any)
├── obsp/                   # Cell-cell graphs (binary sparse)
├── varp/                   # Gene-gene graphs (if any)
├── layers/                 # Additional matrices (binary sparse)
├── raw/                    # Raw data (if present)
│   ├── matrix.data.npy
│   ├── matrix.indices.npy
│   ├── matrix.indptr.npy
│   └── var.parquet
└── uns.json                # Unstructured metadata
```

### v0.1.3 Format Features

- **Binary CSR/CSC format**: Sparse matrices stored as numpy arrays
  - CSR (default): Fastest write from Python (scanpy outputs CSR)
  - CSC: Faster R read performance
  - Format auto-detected from `shape.json` metadata
- **cells×genes orientation**: No transpose needed on Python read
- **Backward compatible**: Reads legacy MTX-based .scio files

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
