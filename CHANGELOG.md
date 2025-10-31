# Changelog

All notable changes to scBridge will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### In Development
- Initial package implementation
- Testing and validation in progress

## [1.0.0] - TBD

### Added
- **Python package (scbridge)**:
  - `write()` function to save AnnData to tar format
  - `read()` function to load AnnData from tar format
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

- **R package (scBridge)**:
  - `write()` function to save Seurat/SCE to tar format
  - `read()` function to load as Seurat/SCE from tar format
  - Full support for Seurat objects
  - Full support for SingleCellExperiment objects
  - Advanced `write_components()` and `read_components()` for direct component manipulation

- **Cross-platform compatibility**:
  - Python (AnnData) ↔ R (Seurat) round-trip preservation
  - Python (AnnData) ↔ R (SingleCellExperiment) round-trip preservation
  - Universal format using MTX and Parquet

- **Optimizations**:
  - Uncompressed tar for fast extraction (~2-3s)
  - Parquet format for metadata (10-50x faster than CSV)
  - MTX format for sparse matrices (universal compatibility)
  - Support for 1M+ cell datasets

- **Documentation**:
  - Complete API documentation for Python and R
  - Format specification (v1.0)
  - Cross-platform usage examples
  - Release guide

### Technical Details
- Python: Requires Python >= 3.8, anndata >= 0.8.0, pandas >= 1.5.0
- R: Requires R >= 4.0.0, Matrix, arrow, jsonlite
- Format: scBridge v1.0 (see specs/format_v1.md)

## Release Notes Template

### [X.Y.Z] - YYYY-MM-DD

#### Added
- New features

#### Changed
- Changes to existing features

#### Deprecated
- Features that will be removed in future versions

#### Removed
- Removed features

#### Fixed
- Bug fixes

#### Security
- Security improvements

---

## Version History (Post-Release)

*Version history will be added after releases*

---

## Planned Features (Future Releases)

### v1.1.0 (Potential)
- [ ] Chunked reading for very large datasets
- [ ] Progress bars for long operations
- [ ] Validation utilities
- [ ] More detailed error messages

### v1.2.0 (Potential)
- [ ] Cloud storage support (S3, GCS)
- [ ] Streaming read for memory-constrained systems
- [ ] Additional metadata preservation

### v2.0.0 (Potential - Breaking Changes)
- [ ] Format v2.0 with enhanced compression options
- [ ] Plugin system for custom components
- [ ] Advanced querying capabilities

---

## Contributing

To contribute to scBridge, please:
1. Check existing issues and pull requests
2. Follow the coding standards
3. Add tests for new features
4. Update documentation
5. Submit a pull request

See [CONTRIBUTING.md](CONTRIBUTING.md) for details.

---

## Questions or Issues?

- GitHub Issues: https://github.com/YOUR_USERNAME/scbridge/issues
- Documentation: https://github.com/YOUR_USERNAME/scbridge
