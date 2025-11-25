#!/usr/bin/env Rscript
# =============================================================================
# R usage examples for scio
#
# Demonstrates how to use scio to save and load Seurat/SingleCellExperiment
# objects and exchange data with Python.
# =============================================================================

# Try loading installed package, otherwise source from local files
if (!requireNamespace("scio", quietly = TRUE)) {
  # Source the package functions for development/testing
  source_dir <- function(path) {
    files <- list.files(path, pattern = "\\.R$", full.names = TRUE)
    for (f in files) source(f)
  }
  source_dir("../R/R")
  message("Note: Using local source files (package not installed)")
} else {
  library(scio)
}

cat("============================================================\n")
cat("scio R Usage Examples\n")
cat("============================================================\n")

# Check which packages are available
has_seurat <- requireNamespace("Seurat", quietly = TRUE)
has_sce <- requireNamespace("SingleCellExperiment", quietly = TRUE)

# =============================================================================
# 1. Load data from .scio format (Python-created)
# =============================================================================

if (has_seurat) {
  cat("\n1. Loading Python-created .scio file as Seurat...\n")
  seurat <- scio_read("sample_data.scio", output = "Seurat")
  cat("   Loaded:", ncol(seurat), "cells x", nrow(seurat), "genes\n")
  cat("   Reductions:", paste(names(seurat@reductions), collapse = ", "), "\n")
} else {
  cat("\n1. Skipping Seurat example (package not installed)\n")
}

if (has_sce) {
  cat("\n2. Loading as SingleCellExperiment...\n")
  sce <- scio_read("sample_data.scio", output = "SCE")
  cat("   Loaded:", ncol(sce), "cells x", nrow(sce), "genes\n")
  cat("   reducedDims:", paste(SingleCellExperiment::reducedDimNames(sce), collapse = ", "), "\n")
} else {
  cat("\n2. Skipping SingleCellExperiment example (package not installed)\n")
}

# =============================================================================
# 3. Save and reload (if Seurat available)
# =============================================================================

if (has_seurat) {
  cat("\n3. Saving Seurat object to .scio format...\n")
  scio_write(seurat, "output_from_r.scio", overwrite = TRUE)
  cat("   Saved: output_from_r.scio\n")

  cat("\n4. Verifying round-trip (R -> scio -> R)...\n")
  seurat_reloaded <- scio_read("output_from_r.scio", output = "Seurat")
  cat("   Original cells:", ncol(seurat), "\n")
  cat("   Reloaded cells:", ncol(seurat_reloaded), "\n")
  if (ncol(seurat) == ncol(seurat_reloaded)) {
    cat("   Round-trip successful!\n")
  }
  unlink("output_from_r.scio", recursive = TRUE)
}

if (has_sce) {
  cat("\n5. Saving SingleCellExperiment to .scio format...\n")
  scio_write(sce, "output_sce.scio", overwrite = TRUE)
  cat("   Saved: output_sce.scio\n")

  cat("\n6. Verifying SCE round-trip...\n")
  sce_reloaded <- scio_read("output_sce.scio", output = "SCE")
  cat("   Original cells:", ncol(sce), "\n")
  cat("   Reloaded cells:", ncol(sce_reloaded), "\n")
  if (ncol(sce) == ncol(sce_reloaded)) {
    cat("   Round-trip successful!\n")
  }
  unlink("output_sce.scio", recursive = TRUE)
}

# =============================================================================
# Cross-platform usage note
# =============================================================================

cat("\n7. Cross-platform compatibility...\n")
cat("   Files saved by R can be read in Python:\n")
cat("   >>> import scio\n")
cat("   >>> adata = scio.read('output_from_r.scio')\n")

cat("\n============================================================\n")
cat("Examples completed successfully!\n")
cat("============================================================\n")
