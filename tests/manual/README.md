# Manual Tests for scio R Package

This directory contains manual test scripts for verifying scio functionality.

## test_write_read.R

Tests the complete write/read round-trip with SingleCellExperiment objects.

**What it tests:**
- Writing SCE objects to .scio format with compression
- Reading .scio files back to SCE objects
- Data integrity (counts, metadata, embeddings)
- No spurious numbered files created
- Proper gzip compression

**To run:**
```bash
cd tests/manual
Rscript test_write_read.R
```

**Expected output:**
- All checks should pass with ✓
- No numbered files (like "3") should be created
- Final message: "=== ✓ ALL TESTS PASSED ==="

## Test Results

This test verifies the fixes for:
1. File descriptor leak causing numbered files
2. MTX file reading/writing with gzip compression
3. TSV file compression handling
4. Complete data round-trip integrity
