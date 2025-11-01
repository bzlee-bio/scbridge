# scBridge - R Package

Cross-platform single-cell RNA-seq data storage for R (Seurat/SingleCellExperiment) and Python (AnnData).

## Installation

### From GitHub

```r
# Install devtools if needed
install.packages("devtools")

# Install scBridge from GitHub
devtools::install_github("bzlee-bio/scbridge", subdir = "R")
```

### Local Development

```r
# From the R directory
devtools::install(".")
```

## Quick Start

```r
library(scBridge)
library(Seurat)

# Save Seurat object to .scb format (works with Python!)
write(seurat_obj, "data.scb")

# Load back as Seurat
seurat_obj <- read("data.scb", output = "Seurat")

# Load as SingleCellExperiment
sce_obj <- read("data.scb", output = "SCE")
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

- **Cross-platform compatibility**: Works seamlessly with Python
  - Load in Python as AnnData objects
  - See `../python/README.md` for Python usage

- **Efficient storage**:
  - MTX format for sparse matrices (universal compatibility)
  - Parquet format for metadata (10-50x faster than CSV)
  - Direct folder access (no extraction needed)

- **Large dataset support**: Optimized for 1M+ cells

## API Reference

### write()

```r
write(object, path, overwrite = FALSE, compress = TRUE)
```

Save Seurat or SingleCellExperiment object to .scb folder.

**Parameters:**
- `object`: Seurat or SingleCellExperiment object
- `path` (character): Output .scb folder path
- `overwrite` (logical): Whether to overwrite existing folder (default: FALSE)
- `compress` (logical): Whether to compress MTX files (default: TRUE)

**Examples:**
```r
library(Seurat)
library(scBridge)

# Save Seurat object
write(seurat_obj, "data.scb")

# Overwrite existing file
write(seurat_obj, "data.scb", overwrite = TRUE)
```

### read()

```r
seurat_obj <- read(path, output = "Seurat")
sce_obj <- read(path, output = "SCE")
```

Load data from .scb folder (or legacy tar archive) as Seurat or SingleCellExperiment.

**Parameters:**
- `path` (character): Path to .scb folder (also supports legacy tar archives)
- `output` (character): Output format - "Seurat" or "SCE" (default: "Seurat")

**Returns:**
- Seurat or SingleCellExperiment object

**Examples:**
```r
# Load as Seurat
seurat_obj <- read("data.scb", output = "Seurat")

# Load as SingleCellExperiment
sce_obj <- read("data.scb", output = "SCE")
```

## File Format

The .scb folder structure:

```
data.scb/
├── manifest.json           # Metadata about saved components
├── matrix.mtx              # Expression matrix (genes × cells)
  ├── barcodes.tsv.gz         # Cell IDs
  ├── features.tsv.gz         # Gene IDs
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
  ├── layers/                 # Additional assays (if any)
  ├── raw/                    # Raw data (if present)
  │   ├── matrix.mtx
  │   ├── features.tsv.gz
  │   └── var.parquet
  └── uns.json                # Unstructured metadata
```

## Cross-Platform Usage

### R to Python

```r
# R: Save from Seurat
library(scBridge)
library(Seurat)

write(seurat_obj, "data.scb")
```

```python
# Python: Load as AnnData
import scbridge as sb

adata = sb.read("data.scb")
```

### Python to R

```python
# Python: Save
import anndata as ad
import scbridge as sb

adata = ad.read_h5ad("data.h5ad")
sb.write(adata, "data.scb")
```

```r
# R: Load as Seurat
library(scBridge)

seurat_obj <- read("data.scb", output = "Seurat")
```

## Requirements

- R >= 4.0.0
- Matrix
- arrow
- jsonlite

**Optional (for object conversion):**
- Seurat (for Seurat objects)
- SingleCellExperiment (for SCE objects)

Install dependencies:

```r
# Required packages
install.packages(c("Matrix", "arrow", "jsonlite"))

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
write_components(components, "data.scb")

# Load components
components <- read_components("data.scb")

# Convert to Seurat
seurat_obj <- create_seurat_from_components(components)
```

## License

MIT License

## Contributing

See [CONTRIBUTING.md](../CONTRIBUTING.md)


