# scBridge Package - File Creation Guide

## Status

✅ Package planning complete
✅ Repository structure created
🔄 Ready to create all package files

## What Will Be Created

### Total Files: ~35 files

Due to conversation length limits, I recommend continuing in a fresh session to create all files.

## Files to Create - Checklist

### Python Package (python/scbridge/) - 12 files
- [ ] `__init__.py` - Package initialization
- [ ] `io.py` - Main API (write, read functions)
- [ ] `writers.py` - Save all 10 components to folder
- [ ] `readers.py` - Load all 10 components from folder
- [ ] `formats.py` - MTX, Parquet, JSON handlers
- [ ] `_version.py` - Version info
- [ ] `../setup.py` - Package setup
- [ ] `../pyproject.toml` - Modern Python packaging
- [ ] `../requirements.txt` - Dependencies
- [ ] `../README.md` - Python package docs
- [ ] `../LICENSE` - MIT License
- [ ] `../.gitignore` - Git ignore

### R Package (R/) - 15 files
- [ ] `DESCRIPTION` - R package metadata
- [ ] `NAMESPACE` - Exports
- [ ] `R/io.R` - Main API (write, read)
- [ ] `R/writers.R` - Save from Seurat/SCE
- [ ] `R/readers.R` - Load to Seurat/SCE
- [ ] `R/seurat.R` - Seurat conversions
- [ ] `R/sce.R` - SCE conversions
- [ ] `R/formats.R` - Format handlers
- [ ] `R/utils.R` - Utilities
- [ ] `man/write.Rd` - Documentation
- [ ] `man/read.Rd` - Documentation
- [ ] `tests/testthat.R` - Test runner
- [ ] `tests/testthat/test-io.R` - Tests
- [ ] `README.md` - R package docs
- [ ] `.Rbuildignore` - R build ignore

### Documentation (specs/, examples/) - 5 files
- [ ] `specs/format_v1.md` - Complete format specification
- [ ] `examples/python_to_r.md` - Python → R workflow
- [ ] `examples/r_to_python.md` - R → Python workflow
- [ ] `CONTRIBUTING.md` - Contribution guide
- [ ] `CODE_OF_CONDUCT.md` - Code of conduct

### Release Files - 3 files
- [ ] `RELEASE_GUIDE.md` - Step-by-step release instructions
- [ ] `CHANGELOG.md` - Version history
- [ ] `.github/workflows/ci.yml` - GitHub Actions CI

## Implementation Strategy

The core implementation reuses your proven code:

### Python Writers (writers.py)
```python
# Adapts from: save_r_compatible_complete.py
def save_to_folder(adata, folder_path):
    # Save X, obs, var (core)
    # Save obsm, varm, obsp, varp (optional)
    # Save layers (optional)
    # Save raw (optional)
    # Save uns (optional)
    # Create manifest.json
```

### Python Readers (readers.py)
```python
def load_from_folder(folder_path):
    # Read manifest.json
    # Load X, obs, var (core)
    # Load all optional components
    # Reconstruct AnnData
    return adata
```

### Python Main API (io.py)
```python
import tarfile
import tempfile

def write(adata, path, overwrite=False):
    # Check overwrite
    # Create temp folder
    # save_to_folder(adata, temp_folder)
    # Create tar from temp folder
    # Clean up temp

def read(path):
    # Create temp folder
    # Extract tar to temp folder
    # load_from_folder(temp_folder)
    # Clean up temp
    # return adata
```

### R Readers (readers.R)
```r
# Adapts from: load_scrnaseq_data_utils.R
load_from_folder <- function(folder_path) {
  # Load all components
  # Return list with X, obs, var, obsm, etc.
}
```

### R Converters (seurat.R, sce.R)
```r
create_seurat_from_components <- function(components) {
  # Convert to Seurat object
}

create_sce_from_components <- function(components) {
  # Convert to SCE object
}
```

### R Main API (io.R)
```r
write <- function(object, path, overwrite = FALSE) {
  # Detect if Seurat or SCE
  # Convert to components
  # Save to temp folder
  # Create tar
  # Clean up
}

read <- function(path, output = "Seurat") {
  # Extract tar
  # Load components
  # Convert to Seurat or SCE
  # Clean up
  # return object
}
```

## Next Steps

**Option 1: Continue in New Session (Recommended)**
- Fresh conversation
- Create all ~35 files systematically
- Test as we go

**Option 2: Create Core Files Now**
- Create minimal working version (10-15 files)
- Test basic functionality
- Complete remaining files later

**Option 3: Provide Implementation Scripts**
- Create Python script to generate all Python files
- Create R script to generate all R files
- Run scripts to create package

## Current Progress

✅ Planning: 100%
✅ Design: 100%
✅ Repository structure: 100%
🔄 File creation: 5% (README, summary docs)
⏳ Implementation: 0%
⏳ Testing: 0%
⏳ Documentation: 20%
⏳ Release prep: 0%

## Recommendation

Start a new session with clear goal:
"Create all scBridge package files based on the plan in /home/bz/bz-nas/pretrained_model/cancer_pretrained/FC_model_organized/scbridge/"

This will allow systematic file creation without conversation length limits.

## Your Existing Code Integration

These files will be adapted:
- `/home/bz/.../save_r_compatible_complete.py` → `python/scbridge/writers.py`
- `/home/bz/.../load_scrnaseq_data_utils.R` → `R/readers.R`

The package adds:
- Tar packaging layer
- Simplified API
- Cross-platform testing
- Documentation
- Release infrastructure

## Ready to Continue?

When ready, in a new session I will:
1. Create all Python package files
2. Create all R package files
3. Create documentation
4. Create tests
5. Create release guide for your approval

Estimated time: 2-3 hours of focused implementation.
