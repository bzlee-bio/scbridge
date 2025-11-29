# scBridge Format Specification v0.1.3

## Overview

scBridge uses a folder-based structure with `.scb` extension. This format is designed for:

- **Cross-platform compatibility**: Works with Python (AnnData) and R (Seurat/SCE)
- **Complete data preservation**: Saves ALL 10 AnnData components
- **Efficient storage**: Uses MTX for sparse matrices, Parquet for DataFrames
- **Direct folder access**: No extraction needed, instant loading
- **Large dataset support**: Optimized for 1M+ cells

## Folder Structure

```
dataset.scb/
├── manifest.json              # Required: Format metadata
├── matrix.mtx                 # Required: Main expression matrix
├── barcodes.tsv[.gz]          # Required: Cell IDs (optionally gzipped)
├── features.tsv[.gz]          # Required: Gene IDs (optionally gzipped)
├── obs.parquet                # Required: Cell metadata
├── var.parquet                # Optional: Gene metadata
├── obsm/                      # Optional: Cell embeddings folder
│   ├── X_pca.parquet
│   ├── X_umap.parquet
│   └── [embedding_name].parquet
├── varm/                      # Optional: Gene embeddings folder
│   └── [embedding_name].parquet
├── obsp/                      # Optional: Cell-cell graphs folder
│   ├── distances.mtx
│   ├── connectivities.mtx
│   └── [graph_name].mtx
├── varp/                      # Optional: Gene-gene graphs folder
│   └── [graph_name].mtx
├── layers/                    # Optional: Additional matrices folder
│   └── [layer_name].mtx
├── raw/                       # Optional: Raw data folder
│   ├── matrix.mtx
│   ├── features.tsv[.gz]      # Optionally gzipped
│   └── var.parquet
└── uns.json                   # Optional: Unstructured metadata
```

**Note:** Files marked with `[.gz]` may be gzipped depending on the `compress` parameter.

## File Formats

### 1. manifest.json (Required)

JSON file containing metadata about the dataset.

**Schema:**

```json
{
  "format": "scBridge v1.0",
  "created_by": "scbridge.write" or "scBridge::write",
  "dimensions": {
    "n_obs": <integer>,  // Number of cells
    "n_vars": <integer>  // Number of genes
  },
  "components": {
    "X": true,
    "obs": true,
    "var": <boolean>,
    "obsm": [<list of embedding names>],
    "varm": [<list of gene embedding names>],
    "obsp": [<list of cell-cell graph names>],
    "varp": [<list of gene-gene graph names>],
    "layers": [<list of layer names>],
    "raw": <boolean>,
    "uns": [<list of uns keys>]
  },
  "files": {
    "X": "matrix.mtx",
    "barcodes": "barcodes.tsv" or "barcodes.tsv.gz",
    "features": "features.tsv" or "features.tsv.gz",
    "obs": "obs.parquet",
    // ... mapping of components to file paths
  }
}
```

### 2. matrix.mtx (Required)

Main expression matrix in Matrix Market format.

**Specifications:**
- Format: Matrix Market (MTX), uncompressed (no gzip for performance)
- Orientation: **genes × cells** (transposed for R compatibility)
- Type: Sparse matrix (CSR/CSC)
- Will be transposed to **cells × genes** when loaded

**Tools:**
- Python: `scipy.io.mmwrite()`, `scipy.io.mmread()`
- R: `Matrix::writeMM()`, `Matrix::readMM()`

**Note:** MTX files are not gzipped by default for performance reasons (gzip is 10-100x slower for large matrices).

### 3. barcodes.tsv[.gz] (Required)

Cell IDs (barcodes) file.

**Specifications:**
- Format: Tab-separated values, optionally gzip compressed
- Columns: Single column with cell IDs
- No header
- Order must match matrix columns
- Extension: `.tsv` or `.tsv.gz` based on `compress` parameter

**Example:**
```
cell_0
cell_1
cell_2
```

### 4. features.tsv[.gz] (Required)

Gene IDs (features) file.

**Specifications:**
- Format: Tab-separated values, optionally gzip compressed
- Columns: 3 columns (gene_id, gene_name, feature_type)
- No header
- Order must match matrix rows
- Extension: `.tsv` or `.tsv.gz` based on `compress` parameter

**Example:**
```
ENSG00001    ENSG00001    Gene Expression
ENSG00002    ENSG00002    Gene Expression
ENSG00003    ENSG00003    Gene Expression
```

### 5. obs.parquet (Required)

Cell metadata DataFrame.

**Specifications:**
- Format: Apache Parquet
- Index: Cell IDs (matching barcodes.tsv.gz)
- Columns: Any cell-level metadata
- Preserves pandas dtypes (Float64, string, etc.)

**Advantages:**
- 10-50x faster than CSV for large data
- Preserves nullable dtypes (Float64, string)
- Column-oriented for efficient queries

### 6. var.parquet (Optional)

Gene metadata DataFrame.

**Specifications:**
- Format: Apache Parquet
- Index: Gene IDs (matching features.tsv.gz)
- Columns: Any gene-level metadata

### 7. obsm/*.parquet (Optional)

Cell embeddings folder.

**Specifications:**
- Each file: One embedding (e.g., X_pca.parquet, X_umap.parquet)
- Format: Apache Parquet
- Index: Cell IDs
- Columns: Named as `<embedding_name>_1`, `<embedding_name>_2`, etc.

**Example obsm/X_pca.parquet:**
```
            X_pca_1    X_pca_2    X_pca_3    ...
cell_0      0.234     -1.567      0.891
cell_1     -0.456      2.123     -0.345
```

### 8. varm/*.parquet (Optional)

Gene embeddings folder.

**Specifications:**
- Each file: One gene embedding
- Format: Apache Parquet
- Index: Gene IDs
- Columns: Named as `<embedding_name>_1`, `<embedding_name>_2`, etc.

### 9. obsp/*.mtx (Optional)

Cell-cell graphs folder.

**Specifications:**
- Each file: One graph (e.g., distances.mtx, connectivities.mtx)
- Format: Matrix Market (MTX), uncompressed
- Orientation: **cells × cells**
- Type: Sparse matrix

**Common graphs:**
- `distances.mtx`: K-nearest neighbor distances
- `connectivities.mtx`: K-nearest neighbor connectivities

### 10. varp/*.mtx (Optional)

Gene-gene graphs folder.

**Specifications:**
- Each file: One gene-gene graph
- Format: Matrix Market (MTX), uncompressed
- Orientation: **genes × genes**
- Type: Sparse matrix

### 11. layers/*.mtx (Optional)

Additional matrices folder.

**Specifications:**
- Each file: One layer/assay (e.g., raw_counts.mtx)
- Format: Matrix Market (MTX), uncompressed
- Orientation: **genes × cells** (transposed for R compatibility)
- Will be transposed to **cells × genes** when loaded

**Common layers:**
- `counts.mtx`: Raw counts
- `logcounts.mtx`: Log-normalized counts

### 12. raw/ (Optional)

Raw data folder (for adata.raw).

**Contents:**
- `matrix.mtx`: Raw expression matrix (genes × cells)
- `features.tsv[.gz]`: Raw gene IDs (optionally gzipped)
- `var.parquet`: Raw gene metadata

**Specifications:**
- Same format as main matrix files
- Typically contains pre-filtered gene set

### 13. uns.json (Optional)

Unstructured metadata.

**Specifications:**
- Format: JSON
- Content: JSON-serializable metadata
- Non-serializable objects are skipped with warning

**Example:**
```json
{
  "experiment_name": "PBMC_10k",
  "sequencing_date": "2024-01-15",
  "n_pcs": 50,
  "neighbors": {
    "n_neighbors": 15,
    "method": "umap"
  }
}
```

## Data Type Mapping

### Python (AnnData) → Storage

| Component | Storage Format | Notes |
|-----------|---------------|-------|
| adata.X | matrix.mtx | Transposed to genes × cells |
| adata.obs | obs.parquet | Preserves pandas dtypes |
| adata.var | var.parquet | Preserves pandas dtypes |
| adata.obsm | obsm/*.parquet | Each key → separate file |
| adata.varm | varm/*.parquet | Each key → separate file |
| adata.obsp | obsp/*.mtx | Each key → separate file |
| adata.varp | varp/*.mtx | Each key → separate file |
| adata.layers | layers/*.mtx | Each key → separate file, transposed |
| adata.raw | raw/ folder | Full raw structure |
| adata.uns | uns.json | JSON-serializable only |

### R (Seurat) → Storage

| Component | Storage Format | Notes |
|-----------|---------------|-------|
| counts (assay) | matrix.mtx | Transposed to genes × cells |
| meta.data | obs.parquet | Metadata as Parquet |
| reductions | obsm/*.parquet | Each reduction → separate file |
| graphs | obsp/*.mtx | Each graph → separate file |
| other assays | layers/*.mtx | Each assay → separate file |
| misc$uns | uns.json | Unstructured metadata |

### R (SingleCellExperiment) → Storage

| Component | Storage Format | Notes |
|-----------|---------------|-------|
| assay(sce, "counts") | matrix.mtx | Transposed to genes × cells |
| colData(sce) | obs.parquet | Cell metadata |
| rowData(sce) | var.parquet | Gene metadata |
| reducedDims(sce) | obsm/*.parquet | Each reducedDim → separate file |
| metadata(sce)$obsp | obsp/*.mtx | Cell-cell graphs |
| metadata(sce)$varm | varm/*.parquet | Gene embeddings |
| metadata(sce)$varp | varp/*.mtx | Gene-gene graphs |
| other assays | layers/*.mtx | Each assay → separate file |
| metadata(sce)$uns | uns.json | Unstructured metadata |

## Compression Strategy

- **Folder format**: Direct access (no archive extraction needed)
- **MTX files**: Uncompressed by default (*.mtx) for performance
- **TSV files**: Optionally gzipped (*.tsv or *.tsv.gz) based on `compress` parameter
- **Parquet files**: Built-in compression (no .gz)
- **JSON files**: Uncompressed (small file size)

**Rationale:**
- Folder format: No extraction overhead, instant access
- Uncompressed MTX: 10-100x faster than gzipped MTX for large matrices
- Optional TSV compression: User can choose between size and speed
- Parquet: Already efficient, built-in compression

## Version Compatibility

**Version 0.1.3 (Beta):**
- Folder-based format (no tar archive)
- All 10 AnnData components supported
- Configurable sparse format: CSR (default, fastest write) or CSC (faster R read)
- Optional gzip compression for TSV files
- Python 3.8+ and R 4.0+ compatible

**Legacy support:**
- Readers support legacy tar archives for backward compatibility
- Writers only create folder format

**Future versions:**
- Backward compatible: Newer readers can read older formats
- Forward compatible: Older readers may not support new features

## Performance Characteristics

**Typical 1M cell dataset:**
- Save time: ~15-20 seconds
- Load time: ~5-8 seconds
- No extraction needed (folder format)
- File size: ~600 MB

**Comparison to alternatives:**
- h5ad: Similar performance, but doesn't preserve all dtypes
- CSV: 10-50x slower for large datasets
- Legacy tar: Added extraction overhead (~2-3s)

## Validation

Readers should validate:

1. **manifest.json exists**: Required file
2. **Required files exist**: matrix.mtx, barcodes.tsv[.gz], features.tsv[.gz], obs.parquet
3. **Dimensions match**: n_obs and n_vars in manifest match actual data
4. **Index alignment**: barcodes match obs index, features match var index

## Error Handling

**Missing required files:**
- Error: "File not found: [filename]"

**Dimension mismatch:**
- Error: "Dimension mismatch: expected [n] cells, found [m]"

**Corrupted files:**
- Error: "Failed to read [filename]: [specific error]"

## Extensions

Future extensions may include:

- Additional metadata in manifest
- Alternative compression algorithms
- Chunked reading for very large datasets
- Cloud-optimized formats

## References

- **Matrix Market format**: https://math.nist.gov/MatrixMarket/formats.html
- **Apache Parquet**: https://parquet.apache.org/
- **AnnData**: https://anndata.readthedocs.io/
