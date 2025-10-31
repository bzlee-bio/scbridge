# scBridge Release Guide

## IMPORTANT: DO NOT RELEASE WITHOUT USER APPROVAL

This guide is for **REFERENCE ONLY**. Do NOT execute any release steps without explicit approval from the repository owner.

## Pre-Release Checklist

### 1. Testing

**Python Package:**
```bash
cd scbridge/python

# Create test environment
python -m venv test_env
source test_env/bin/activate

# Install package
pip install -e .

# Test basic functionality
python -c "
import anndata as ad
import scbridge as sb
import numpy as np

# Create test data
adata = ad.AnnData(np.random.rand(100, 50))
adata.obs['cell_type'] = ['A', 'B'] * 50
adata.obsm['X_pca'] = np.random.rand(100, 10)

# Test write
sb.write(adata, 'test.tar')
print('Write: OK')

# Test read
adata2 = sb.read('test.tar')
print(f'Read: OK - Shape {adata2.shape}')
"
```

**R Package:**
```r
# In R console
cd scbridge/R

# Install dependencies
install.packages(c("Matrix", "arrow", "jsonlite", "devtools"))

# Install package
devtools::install(".")

# Test basic functionality
library(scBridge)
library(Matrix)

# Create test Seurat-like components
X <- Matrix::rsparsematrix(100, 50, density = 0.3)
rownames(X) <- paste0("cell_", 1:100)
colnames(X) <- paste0("gene_", 1:50)

obs <- data.frame(
  cell_type = rep(c("A", "B"), 50),
  row.names = rownames(X)
)

var <- data.frame(row.names = colnames(X))

components <- list(
  X = X,
  obs = obs,
  var = var,
  obsm = list(),
  varm = list(),
  obsp = list(),
  varp = list(),
  layers = list(),
  uns = list()
)

# Test write
write_components(components, "test.tar")
print("Write: OK")

# Test read
components2 <- read_components("test.tar")
print(paste("Read: OK - Shape", nrow(components2$X), "x", ncol(components2$X)))
```

### 2. Cross-Platform Testing

**Python → R:**
```python
# Python: Save
import anndata as ad
import scbridge as sb
adata = ad.read_h5ad("real_data.h5ad")
sb.write(adata, "test_python_to_r.tar")
```

```r
# R: Load
library(scBridge)
seurat_obj <- read("test_python_to_r.tar", output = "Seurat")
print(seurat_obj)
```

**R → Python:**
```r
# R: Save
library(scBridge)
library(Seurat)
write(seurat_obj, "test_r_to_python.tar")
```

```python
# Python: Load
import scbridge as sb
adata = sb.read("test_r_to_python.tar")
print(adata)
```

### 3. Documentation Review

- [ ] README.md files are up-to-date
- [ ] API documentation is complete
- [ ] Examples work correctly
- [ ] CHANGELOG.md is updated
- [ ] Version numbers are consistent

### 4. Code Quality

**Python:**
```bash
# Format check (optional)
cd scbridge/python
python -m pylint scbridge/

# Type checking (optional)
python -m mypy scbridge/
```

**R:**
```r
# Check package
devtools::check()

# Build documentation
devtools::document()
```

## Release Steps (REFERENCE ONLY - DO NOT EXECUTE)

### Python Package Release

**1. Update Version:**
```bash
# Edit scbridge/python/scbridge/_version.py
__version__ = "1.0.0"

# Edit scbridge/python/pyproject.toml
version = "1.0.0"
```

**2. Build Package:**
```bash
cd scbridge/python

# Install build tools
pip install build twine

# Build distributions
python -m build

# This creates:
# - dist/scbridge-1.0.0.tar.gz
# - dist/scbridge-1.0.0-py3-none-any.whl
```

**3. Test Upload (TestPyPI):**
```bash
# Upload to TestPyPI first
python -m twine upload --repository testpypi dist/*

# Test installation
pip install --index-url https://test.pypi.org/simple/ scbridge
```

**4. Production Upload (PyPI):**
```bash
# ONLY after explicit approval
python -m twine upload dist/*

# Package will be available at:
# https://pypi.org/project/scbridge/
```

**5. GitHub Release:**
```bash
# Tag the release
git tag -a python-v1.0.0 -m "Python package v1.0.0"
git push origin python-v1.0.0

# Create GitHub release with:
# - Release notes from CHANGELOG.md
# - Attach dist/*.tar.gz and dist/*.whl files
```

### R Package Release

**1. Update Version:**
```r
# Edit scbridge/R/DESCRIPTION
Version: 1.0.0
```

**2. Build Package:**
```r
cd scbridge/R

# Build package
devtools::build()

# This creates:
# - scBridge_1.0.0.tar.gz
```

**3. Check Package:**
```r
# Run R CMD check
devtools::check()

# Should pass with 0 errors, 0 warnings
```

**4. GitHub Release:**
```bash
# Tag the release
git tag -a r-v1.0.0 -m "R package v1.0.0"
git push origin r-v1.0.0

# Create GitHub release
```

**5. CRAN Submission (Optional - requires user decision):**

CRAN submission is a lengthy process. Only proceed if explicitly requested.

```r
# Build for CRAN
devtools::check(cran = TRUE)

# Submit to CRAN (ONLY after approval)
devtools::submit_cran()
```

### Documentation Website (Optional)

If creating a documentation website:

```bash
# Python: Use Sphinx or MkDocs
# R: Use pkgdown

cd scbridge/R
Rscript -e "pkgdown::build_site()"

# Deploy to GitHub Pages
```

## Post-Release

### 1. Announcement

Prepare announcement template:

```markdown
# scBridge v1.0.0 Released

We're excited to announce the release of scBridge v1.0.0!

## Features
- Cross-platform single-cell data storage (Python ↔ R)
- Complete preservation of all AnnData components
- Fast and efficient for large datasets (1M+ cells)

## Installation

Python:
pip install scbridge

R:
devtools::install_github("USERNAME/scbridge", subdir = "R")

## Links
- GitHub: https://github.com/USERNAME/scbridge
- Documentation: [link]
- PyPI: https://pypi.org/project/scbridge/

## Feedback
Please report issues at: https://github.com/USERNAME/scbridge/issues
```

### 2. Update README

Update main README.md with:
- Installation instructions
- Links to released versions
- Citation information

### 3. Monitor Issues

- Watch for bug reports
- Respond to user questions
- Plan for patch releases if needed

## Version Numbering

Follow Semantic Versioning (semver):

- **MAJOR.MINOR.PATCH** (e.g., 1.0.0)
- **MAJOR**: Breaking changes
- **MINOR**: New features (backward compatible)
- **PATCH**: Bug fixes

Examples:
- 1.0.0 → 1.0.1: Bug fix
- 1.0.0 → 1.1.0: New feature added
- 1.0.0 → 2.0.0: Breaking API change

## Emergency Rollback

If a critical bug is found after release:

**Python:**
```bash
# Yank the release from PyPI (makes it unavailable for new installs)
python -m twine upload --repository pypi --skip-existing --skip-existing scbridge-1.0.0*

# Release fixed version immediately as 1.0.1
```

**R:**
```r
# Update GitHub release with warning
# Release fixed version as 1.0.1
```

## Checklist Summary

Before ANY release:

- [ ] All tests pass
- [ ] Cross-platform compatibility verified
- [ ] Documentation is complete and accurate
- [ ] CHANGELOG.md is updated
- [ ] Version numbers are updated consistently
- [ ] User has explicitly approved release
- [ ] GitHub repository is set up
- [ ] PyPI/CRAN accounts are ready (if needed)

## User Approval Required

Before executing ANY of the above steps:

1. Present complete release plan to user
2. Get explicit written approval
3. Confirm GitHub repository URL
4. Confirm PyPI package name
5. Confirm any CRAN submission plans

## Support Contact

For questions about this release guide:
- Review this guide with the repository owner
- Do not proceed without explicit approval
- Discuss any uncertainties before taking action
