# Changelog

All notable changes to scio will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.3] - 2025-01-29

### Added
- **Configurable sparse format**: New `sparse_format` parameter for write functions
  - CSR (default): Fastest write from Python (scanpy outputs CSR)
  - CSC: Faster R read performance
  - Readers auto-detect format from `shape.json` metadata

### Changed
- Default sparse format is now CSR (was CSC in v0.1.2)
- Python `scio.write()` accepts `sparse_format='csr'` or `sparse_format='csc'`
- R `scio_write()` accepts `sparse_format = "csr"` or `sparse_format = "csc"`

### API Changes
- Python: `scio.write(adata, path, sparse_format='csr', overwrite=False, update=False)`
- R: `scio_write(object, path, sparse_format = "csr", overwrite = FALSE, update = FALSE)`

## [0.1.2] - 2025-01-25

### Added
- **Binary CSC format**: New high-performance sparse matrix storage
  - Stores matrices as binary numpy files (data.npy, indices.npy, indptr.npy)
  - 3.4x faster than H5AD in R, 14.2x faster than MTX
  - Zero-copy transpose using MatrixExtra::t_shallow() in R

- **Orientation metadata**: `cells_x_genes` orientation for v0.1.2+ format
  - Expression matrices stored in cells x genes orientation
  - Eliminates transpose on Python read
  - R readers check orientation and transpose only when needed

- **reticulate+numpy loader**: Primary .npy file loader in R
  - Replaces RcppCNPy which had segfault issues on large matrices
  - Falls back to RcppCNPy if reticulate/numpy unavailable

### Changed
- Expression matrices (X, layers, raw.X) now stored in cells x genes orientation
- Manifest includes `orientation` field to indicate matrix layout
- R readers detect format version and handle transpose accordingly

### Performance
- Python read: ~8.9s (vs H5AD 6.1s)
- Python write: ~7.8s (vs H5AD 9.5s)
- R read: ~29.4s (vs H5AD/zellkonverter 99.3s, MTX 418.4s)

## [0.1.1] - 2025-01-20

### Added
- **Python package (scio)**:
  - `write()` function to save AnnData to .scio folder format
  - `read()` function to load AnnData from .scio folder
  - Complete preservation of all 10 AnnData components:
    - X (expression matrix)
    - obs (cell metadata)
    - var (gene metadata)
    - obsm (cell embeddings)
    - varm (gene embeddings)
    - obsp (cell-cell graphs)
    - varp (gene-gene graphs)
    - layers (additional matrices)
    - raw (raw counts)
    - uns (unstructured metadata)
  - Incremental updates with hash-based change detection

- **R package (scio)**:
  - `scio_write()` function to save Seurat/SCE to .scio folder
  - `scio_read()` function to load as Seurat/SCE
  - Full support for Seurat objects
  - Full support for SingleCellExperiment objects

- **Cross-platform compatibility**:
  - Python (AnnData) <-> R (Seurat) round-trip preservation
  - Python (AnnData) <-> R (SingleCellExperiment) round-trip preservation
  - Universal format using MTX and Parquet

- **Efficient storage**:
  - Folder format (no tar extraction needed)
  - Parquet format for metadata (10-50x faster than CSV)
  - MTX format for sparse matrices
  - Optional gzip compression for MTX files

### Technical Details
- Python: Requires Python >= 3.8, anndata >= 0.8.0, pandas >= 1.5.0
- R: Requires R >= 4.0.0, Matrix, arrow, jsonlite

---

## Planned Features (Future Releases)

- [ ] Chunked reading for very large datasets
- [ ] Progress bars for long operations
- [ ] Cloud storage support (S3, GCS)
- [ ] Streaming read for memory-constrained systems

---

## Contributing

To contribute to scio, please:
1. Check existing issues and pull requests
2. Follow the coding standards
3. Add tests for new features
4. Update documentation
5. Submit a pull request

See [CONTRIBUTING.md](CONTRIBUTING.md) for details.

---

## Questions or Issues?

- GitHub Issues: https://github.com/bzlee-bio/scio/issues
- Documentation: https://github.com/bzlee-bio/scio
