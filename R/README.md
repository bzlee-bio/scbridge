# scio - R Package (v0.1.3)

Cross-platform single-cell RNA-seq data storage for R (Seurat/SingleCellExperiment) and Python (AnnData).

**v0.1.3 Highlights:**
- Binary CSR/CSC format for fast sparse matrix I/O
- Configurable sparse format: CSR (fastest write) or CSC (faster R read)
- 3.4x faster than H5AD (zellkonverter), 14.2x faster than MTX
- Zero-copy transpose using MatrixExtra::t_shallow()

## Installation

### From GitHub

```r
# Install devtools if needed
install.packages("devtools")

# Install scio from GitHub
devtools::install_github("bzlee-bio/scio", subdir = "R")
```

### Local Development

```r
# From the R directory
devtools::install(".")
```

## Quick Start

```r
library(scio)
library(Seurat)

# Save Seurat object to .scio format (works with Python!)
scio_write(seurat_obj, "data.scio")

# Or use CSC format for faster R reading
scio_write(seurat_obj, "data.scio", sparse_format = "csc")

# Load back as Seurat
seurat_obj <- scio_read("data.scio", output = "Seurat")

# Load as SingleCellExperiment
sce_obj <- scio_read("data.scio", output = "SCE")
```

## Features

- **Complete data preservation**: Saves ALL components
  - Expression matrices
  - Cell metadata (colData)
  - Gene metadata (rowData)
  - Embeddings (reducedDims/reductions)
  - Cell-cell graphs
  - Gene embeddings
  - Layers/assays
  - Raw counts
  - Unstructured metadata

- **Fast binary CSR/CSC format** (v0.1.3):
  - Sparse matrices stored as numpy arrays (.npy files)
  - CSR (default): Fastest write from Python (scanpy outputs CSR)
  - CSC: Faster R read performance
  - Uses reticulate+numpy for reading .npy files
  - Zero-copy transpose with MatrixExtra::t_shallow()
  - Benchmark: 29.4s vs H5AD 99.3s vs MTX 418.4s

- **Incremental updates**: Hash-based change detection
  - Only rewrites modified components
  - Dramatically faster for iterative workflows
  - Use `update = TRUE` parameter

- **Cross-platform compatibility**: Works seamlessly with Python
  - Load in Python as AnnData objects
  - See `../python/README.md` for Python usage

- **Efficient storage**:
  - Binary numpy arrays for sparse matrices
  - Parquet format for metadata
  - Direct folder access (no extraction needed)

- **Large dataset support**: Optimized for 1M+ cells

## API Reference

### scio_write()

```r
scio_write(object, path, sparse_format = "csr", overwrite = FALSE, update = FALSE)
```

Save Seurat or SingleCellExperiment object to .scio folder using binary sparse format.
Data is stored in cells × genes orientation for fast cross-platform loading.

**Parameters:**
- `object`: Seurat or SingleCellExperiment object
- `path` (character): Output .scio folder path
- `sparse_format` (character): Format for sparse matrices - `"csr"` (default, fastest write) or `"csc"` (faster R read)
- `overwrite` (logical): Whether to overwrite existing folder (default: FALSE)
- `update` (logical): Whether to perform incremental update (default: FALSE). Only changed components will be rewritten.

**Examples:**
```r
library(Seurat)
library(scio)

# Save Seurat object
scio_write(seurat_obj, "data.scio")

# Overwrite existing file
scio_write(seurat_obj, "data.scio", overwrite = TRUE)

# Use CSC format for faster R reading
scio_write(seurat_obj, "data.scio", sparse_format = "csc")

# Incremental update (only writes changed components)
scio_write(seurat_obj, "data.scio", update = TRUE)
```

### scio_read()

```r
seurat_obj <- scio_read(path, output = "Seurat")
sce_obj <- scio_read(path, output = "SCE")
```

Load data from .scio folder (or legacy tar archive) as Seurat or SingleCellExperiment.

**Parameters:**
- `path` (character): Path to .scio folder (also supports legacy tar archives)
- `output` (character): Output format - "Seurat" or "SCE" (default: "Seurat")

**Returns:**
- Seurat or SingleCellExperiment object

**Examples:**
```r
# Load as Seurat
seurat_obj <- scio_read("data.scio", output = "Seurat")

# Load as SingleCellExperiment
sce_obj <- scio_read("data.scio", output = "SCE")
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
├── layers/                 # Additional assays (binary sparse)
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
- **cells×genes orientation**: R uses MatrixExtra::t_shallow() for zero-copy transpose
- **Backward compatible**: Reads legacy MTX-based .scio files

## Cross-Platform Usage

### R to Python

```r
# R: Save from Seurat
library(scio)
library(Seurat)

scio_write(seurat_obj, "data.scio")
```

```python
# Python: Load as AnnData
import scio

adata = scio.read("data.scio")
```

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

## Requirements

- R >= 4.0.0
- Matrix
- MatrixExtra (for efficient transpose)
- RcppCNPy (for writing .npy files)
- arrow
- jsonlite
- digest (for incremental updates)
- reticulate + numpy (for reading .npy files)

**Optional (for object conversion):**
- Seurat (for Seurat objects)
- SingleCellExperiment (for SCE objects)

Install dependencies:

```r
# Required packages
install.packages(c("Matrix", "MatrixExtra", "RcppCNPy", "arrow", "jsonlite", "digest", "reticulate"))

# Set up Python environment with numpy
reticulate::py_install("numpy")

# Optional packages
install.packages("Seurat")

# For SingleCellExperiment (via Bioconductor)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
```

## Advanced Usage

### Working with Components Directly

```r
# Extract components from Seurat object
components <- extract_components_from_seurat(seurat_obj)

# Modify components
components$uns$custom_info <- "my_data"

# Save components
write_components(components, "data.scio")

# Load components
components <- read_components("data.scio")

# Convert to Seurat
seurat_obj <- create_seurat_from_components(components)
```

## License

MIT License

## Contributing

See [CONTRIBUTING.md](../CONTRIBUTING.md)


