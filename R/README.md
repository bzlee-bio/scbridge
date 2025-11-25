# scio - R Package

Cross-platform single-cell RNA-seq data storage for R (Seurat/SingleCellExperiment) and Python (AnnData).

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

- **Incremental updates**: Hash-based change detection
  - Only rewrites modified components
  - Dramatically faster for iterative workflows
  - Use `update = TRUE` parameter

- **Cross-platform compatibility**: Works seamlessly with Python
  - Load in Python as AnnData objects
  - See `../python/README.md` for Python usage

- **Efficient storage**:
  - MTX format for sparse matrices (universal compatibility)
  - Parquet format for metadata (10-50x faster than CSV)
  - Direct folder access (no extraction needed)

- **Large dataset support**: Optimized for 1M+ cells

## API Reference

### scio_write()

```r
scio_write(object, path, overwrite = FALSE, update = FALSE, compress = TRUE)
```

Save Seurat or SingleCellExperiment object to .scio folder.

**Parameters:**
- `object`: Seurat or SingleCellExperiment object
- `path` (character): Output .scio folder path
- `overwrite` (logical): Whether to overwrite existing folder (default: FALSE)
- `update` (logical): Whether to perform incremental update (default: FALSE). Only changed components will be rewritten.
- `compress` (logical): Whether to gzip compress MTX matrix files (default: TRUE). Compression reduces file size by ~3-5x but is slower to write. TSV files are always gzipped.

**Examples:**
```r
library(Seurat)
library(scio)

# Save Seurat object
scio_write(seurat_obj, "data.scio")

# Overwrite existing file
scio_write(seurat_obj, "data.scio", overwrite = TRUE)

# Incremental update (only writes changed components)
scio_write(seurat_obj, "data.scio", update = TRUE)

# Disable MTX compression for faster writes (larger files)
scio_write(seurat_obj, "data.scio", compress = FALSE)
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

The .scio folder structure:

```
data.scio/
├── manifest.json           # Metadata about saved components
├── matrix.mtx[.gz]         # Expression matrix (optionally gzipped)
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
├── layers/                 # Additional assays (if any)
├── raw/                    # Raw data (if present)
│   ├── matrix.mtx[.gz]
│   ├── features.tsv.gz
│   └── var.parquet
└── uns.json                # Unstructured metadata
```

### Compression

The `compress` parameter controls **MTX file compression only**:

| File Type | compress=TRUE | compress=FALSE |
|-----------|---------------|----------------|
| MTX (matrix, graphs, layers) | `.mtx.gz` (smaller, slower) | `.mtx` (larger, faster) |
| TSV (barcodes, features) | `.tsv.gz` (always) | `.tsv.gz` (always) |
| Parquet (metadata, embeddings) | Built-in compression | Built-in compression |
| JSON (uns) | No compression | No compression |

**When to use `compress=FALSE`:**
- Large datasets where write speed is critical
- Temporary files that will be deleted soon
- When disk space is not a concern

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
- arrow
- jsonlite
- digest (for incremental updates)

**Optional (for object conversion):**
- Seurat (for Seurat objects)
- SingleCellExperiment (for SCE objects)

Install dependencies:

```r
# Required packages
install.packages(c("Matrix", "arrow", "jsonlite", "digest"))

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


