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
import scio as sb

# Load your data
adata = ad.read_h5ad("data.h5ad")

# Save to .scio format (works with R!)
sb.write(adata, "data.scio")

# Load back
adata = sb.read("data.scio")
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
sb.write(adata, path, overwrite=False, compress=True)
```

Save AnnData object to .scio folder.

**Parameters:**
- `adata` (AnnData): AnnData object to save
- `path` (str or Path): Output .scio folder path
- `overwrite` (bool): Whether to overwrite existing folder (default: False)
- `compress` (bool): Whether to compress MTX files (default: True)

**Example:**
```python
sb.write(adata, "data.scio")
sb.write(adata, "data.scio", overwrite=True)
```

### read()

```python
adata = sb.read(path)
```

Load AnnData object from .scio folder (or legacy tar archive).

**Parameters:**
- `path` (str or Path): Path to .scio folder (also supports legacy tar archives)

**Returns:**
- `AnnData`: Reconstructed AnnData object

**Example:**
```python
adata = sb.read("data.scio")
```

## File Format

The .scio folder structure:

```
data.scio/
├── manifest.json           # Metadata about saved components
├── matrix.mtx              # X (expression matrix, genes × cells)
├── barcodes.tsv[.gz]       # Cell IDs (optionally gzipped)
├── features.tsv[.gz]       # Gene IDs (optionally gzipped)
├── obs.parquet             # Cell metadata
├── var.parquet             # Gene metadata
├── obsm/                   # Cell embeddings
│   ├── X_pca.parquet
│   └── X_umap.parquet
├── varm/                   # Gene embeddings (if any)
├── obsp/                   # Cell-cell graphs
│   ├── distances.mtx
│   └── connectivities.mtx
├── varp/                   # Gene-gene graphs (if any)
├── layers/                 # Additional matrices (if any)
├── raw/                    # Raw data (if present)
│   ├── matrix.mtx
│   ├── features.tsv[.gz]   # Optionally gzipped
│   └── var.parquet
└── uns.json                # Unstructured metadata
```

**Note:** Files with `[.gz]` may be gzipped depending on the `compress` parameter used during save.

## Cross-Platform Usage

### Python to R

```python
# Python: Save
import anndata as ad
import scio as sb

adata = ad.read_h5ad("data.h5ad")
sb.write(adata, "data.scio")
```

```r
# R: Load as Seurat
library(scBridge)

seurat_obj <- read("data.scio", output = "Seurat")
```

### R to Python

```r
# R: Save from Seurat
library(scBridge)
library(Seurat)

write(seurat_obj, "data.scio")
```

```python
# Python: Load
import scio as sb

adata = sb.read("data.scio")
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
