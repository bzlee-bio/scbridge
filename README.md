# scio

**Single-cell data I/O for Python and R**

[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

scio provides a universal `.scio` folder format for storing single-cell RNA-seq data that works seamlessly between Python (AnnData/scanpy) and R (Seurat/SingleCellExperiment).

**Key Features:**
- ✅ **Complete data preservation** - All 10 AnnData components (X, obs, var, obsm, varm, obsp, varp, layers, uns, raw)
- ✅ **Cross-platform** - Python ↔ R seamless conversion
- ✅ **Fast** - Binary CSR/CSC format, 3.4x faster than H5AD in R, 14.2x faster than MTX
- ✅ **Efficient** - Binary numpy arrays + Parquet format
- ✅ **Incremental updates** - Hash-based change detection, only rewrites modified components
- ✅ **Simple API** - Just `scio_write()` and `scio_read()` in R, `scio.write()` and `scio.read()` in Python
- ✅ **Easy to inspect** - Standard folder structure with readable formats

## Installation

### Python
```bash
pip install git+https://github.com/bzlee-bio/scio.git#subdirectory=python
```

### R
```R
# Install devtools if needed
install.packages("devtools")

# Install scio from GitHub
devtools::install_github("bzlee-bio/scio", subdir = "R")
```

## Quick Start

### Python → .scio → R

**Python (save):**
```python
import scanpy as sc
import scio

# Load and process data
adata = sc.read_h5ad("pbmc.h5ad")
adata.raw = adata.copy()  # Store raw counts
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Save to .scio (includes ALL components)
scio.write(adata, "pbmc_data.scio")

# Or use CSC format for faster R reading
scio.write(adata, "pbmc_data.scio", sparse_format='csc')
```

**R (load):**
```R
library(scio)
library(Seurat)

# Load as Seurat
seurat <- scio_read("pbmc_data.scio", output = "Seurat")
print(seurat)

# All components preserved:
# - Expression data
# - Raw counts
# - PCA, UMAP embeddings
# - Neighbor graphs
```

### R → .scio → Python

**R (save):**
```R
library(scio)
library(Seurat)

# Process data
seurat <- CreateSeuratObject(counts = pbmc.data)
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)
seurat <- RunUMAP(seurat, dims = 1:30)

# Save to .scio
scio_write(seurat, "pbmc_data.scio")

# Or use CSC format for faster R reading
scio_write(seurat, "pbmc_data.scio", sparse_format = "csc")
```

**Python (load):**
```python
import scio

# Load as AnnData
adata = scio.read("pbmc_data.scio")
print(adata)

# All Seurat reductions loaded to obsm
# All metadata preserved
```

## API Reference

### Python

#### `scio.write(adata, path, sparse_format='csr', overwrite=False, update=False)`
Save AnnData to .scio folder using binary sparse format.

**Parameters:**
- `adata` (AnnData): AnnData object to save
- `path` (str): Output .scio folder path
- `sparse_format` (str): Format for sparse matrices - `'csr'` (default, fastest write) or `'csc'` (faster R read)
- `overwrite` (bool): Whether to overwrite existing folder (default: False)
- `update` (bool): Whether to perform incremental update using hash-based change detection (default: False). Only changed components will be rewritten.

**Example:**
```python
import scio
scio.write(adata, "data.scio")
scio.write(adata, "data.scio", overwrite=True)

# Use CSC format for faster R reading
scio.write(adata, "data.scio", sparse_format='csc')

# Incremental update (only writes changed components)
scio.write(adata, "data.scio", update=True)
```

#### `scio.read(path)`
Load AnnData from .scio folder (or legacy tar archive).

**Parameters:**
- `path` (str): Path to .scio folder (also supports legacy tar archives for backward compatibility)

**Returns:**
- `adata` (AnnData): Loaded AnnData object

**Example:**
```python
adata = scio.read("data.scio")
```

### R

#### `scio_write(object, path, sparse_format = "csr", overwrite = FALSE, update = FALSE)`
Save Seurat or SingleCellExperiment to .scio folder using binary sparse format.

**Parameters:**
- `object`: Seurat or SingleCellExperiment object
- `path` (character): Output .scio folder path
- `sparse_format` (character): Format for sparse matrices - `"csr"` (default, fastest write) or `"csc"` (faster R read)
- `overwrite` (logical): Whether to overwrite existing folder
- `update` (logical): Whether to perform incremental update (default: FALSE). Only changed components will be rewritten.

**Example:**
```R
scio_write(seurat, "data.scio")
scio_write(sce, "data.scio", overwrite = TRUE)

# Use CSC format for faster R reading
scio_write(seurat, "data.scio", sparse_format = "csc")

# Incremental update (only writes changed components)
scio_write(seurat, "data.scio", update = TRUE)
```

#### `scio_read(path, output = c("Seurat", "SCE"))`
Load data from .scio folder.

**Parameters:**
- `path` (character): Path to .scio folder
- `output` (character): Output format ("Seurat" or "SCE")

**Returns:**
- Seurat or SingleCellExperiment object

**Example:**
```R
seurat <- scio_read("data.scio", output = "Seurat")
sce <- scio_read("data.scio", output = "SCE")
```

## What Gets Saved?

scio preserves **ALL** AnnData components:

| Component | Description | Python | R (Seurat) | R (SCE) |
|-----------|-------------|--------|------------|---------|
| X | Expression matrix | adata.X | @assays$RNA@counts | assay(sce, "counts") |
| obs | Cell metadata | adata.obs | @meta.data | colData(sce) |
| var | Gene metadata | adata.var | @assays$RNA@meta.features | rowData(sce) |
| obsm | Cell embeddings | adata.obsm | @reductions | reducedDims(sce) |
| obsp | Cell graphs | adata.obsp | @graphs | - |
| layers | Additional matrices | adata.layers | @assays$RNA@data | assays(sce) |
| uns | Unstructured metadata | adata.uns | @misc | metadata(sce) |
| raw | Raw counts | adata.raw | @assays$raw | altExp(sce) |
| varm | Gene embeddings | adata.varm | @misc | metadata(sce) |
| varp | Gene graphs | adata.varp | @misc | metadata(sce) |

## File Format

scio v0.1.3 uses a `.scio` folder with binary CSR/CSC sparse format:

```
data.scio/
├── manifest.json           # Metadata (format version, orientation, shape)
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
├── obsp/                   # Cell-cell graphs (binary sparse)
│   ├── connectivities.data.npy
│   ├── connectivities.indices.npy
│   └── connectivities.indptr.npy
├── varm/                   # Gene embeddings (if any)
├── varp/                   # Gene-gene graphs (if any)
├── layers/                 # Additional matrices (binary sparse)
├── raw/                    # Raw counts (if any)
│   ├── matrix.data.npy
│   ├── matrix.indices.npy
│   ├── matrix.indptr.npy
│   └── var.parquet
└── uns.json                # Unstructured metadata
```

### v0.1.3 Format Features

- **Binary CSR/CSC format**: Sparse matrices stored as numpy arrays (.npy files)
  - CSR (default): Fastest write from Python (scanpy outputs CSR)
  - CSC: Faster R read performance
  - Format auto-detected from `shape.json` metadata
- **cells×genes orientation**: Expression matrices stored in cells×genes layout
  - No transpose needed on Python read
  - R uses efficient `MatrixExtra::t_shallow()` for zero-copy transpose
- **Backward compatible**: Reads legacy MTX-based .scio files

**Formats used:**
- **Binary NPY** - Sparse matrices (fast binary numpy arrays)
- **Parquet** - Tabular data (metadata, embeddings)
- **JSON** - Unstructured metadata, shape info
- **TSV.gz** - Cell/gene IDs

**Benefits:**
- ✅ 3.4x faster than H5AD in R
- ✅ 14.2x faster than MTX in R
- ✅ Direct memory mapping possible
- ✅ Easy to inspect and debug

## Performance

Benchmark results for 100K cells × 36K genes dataset (v0.1.3):

### Python
| Operation | scio | H5AD |
|-----------|------|------|
| Write | **7.8s** | 9.5s |
| Read | 8.9s | 6.1s |

### R (reading Python-saved data)
| Format | Read Time | vs scio |
|--------|-----------|---------|
| **scio v0.1.3** | **29.4s** | 1x |
| H5AD (zellkonverter) | 99.3s | 3.4x slower |
| MTX | 418.4s | 14.2x slower |

scio v0.1.3 uses binary CSR/CSC format with numpy arrays for maximum performance.

## Requirements

### Python
- Python ≥ 3.8
- anndata ≥ 0.9.0
- pandas ≥ 1.5.0
- numpy ≥ 1.23.0
- scipy ≥ 1.9.0
- pyarrow ≥ 10.0.0

### R
- R ≥ 4.0
- Matrix
- MatrixExtra (for efficient transpose)
- RcppCNPy (for writing .npy files)
- arrow
- jsonlite
- reticulate + numpy (for reading .npy files)
- Seurat (>= 4.0) or SingleCellExperiment (optional)

## Backward Compatibility

**Old `.scb` files from scBridge are fully compatible!**

If you have existing `.scb` files created with the old `scBridge` package, they work seamlessly with `scio`:

```r
# Old .scb files work without any changes
sce <- scio_read("old_data.scb", output = "SCE")

# Or rename for clarity (optional)
# mv old_data.scb old_data.scio
```

The internal format is identical - only the package name and function names changed.

## Examples

See [examples/](examples/) directory for:
- Python → R workflows
- R → Python workflows
- Large dataset handling
- Advanced usage
<!--
## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md). -->

## License

MIT License - see [LICENSE](LICENSE)
<!-- 
## Citation

If you use scio in your research, please cite:

```
@software{scio2025,
  author = {Your Name},
  title = {scio: Cross-platform single-cell RNA-seq data storage},
  year = {2025},
  url = {https://github.com/yourusername/scio}
}
``` -->
<!-- 
## Support

- **Issues**: [GitHub Issues](https://github.com/yourusername/scio/issues)
- **Documentation**: [Read the Docs](https://scio.readthedocs.io)
- **Questions**: [Discussions](https://github.com/yourusername/scio/discussions) -->
