#!/usr/bin/env Rscript
# Manual test for scio with Seurat objects
# Tests that scio_write() and scio_read() work correctly with Seurat

library(Seurat)
library(Matrix)

# Source the package functions directly (for development testing)
source_dir <- function(path) {
  files <- list.files(path, pattern = "\\.R$", full.names = TRUE)
  for (f in files) {
    source(f)
  }
}

# Source all R files from the package
pkg_path <- "../../R/R"
message("Sourcing R files from: ", pkg_path)
source_dir(pkg_path)

# Create a small test Seurat object
set.seed(42)
counts <- Matrix(round(runif(200, 0, 10)), nrow = 20, ncol = 10, sparse = TRUE)
rownames(counts) <- paste0("Gene", 1:20)
colnames(counts) <- paste0("Cell", 1:10)

# Create Seurat object
seurat <- CreateSeuratObject(
  counts = counts,
  project = "TestProject",
  min.cells = 0,
  min.features = 0
)

# Add metadata
seurat$sample <- rep(c("A", "B"), each = 5)
seurat$celltype <- rep(c("T", "B", "NK"), length.out = 10)

# Add PCA reduction
pca_embeddings <- matrix(rnorm(10 * 5), ncol = 5)
rownames(pca_embeddings) <- colnames(seurat)
colnames(pca_embeddings) <- paste0("PC_", 1:5)
seurat[["pca"]] <- CreateDimReducObject(
  embeddings = pca_embeddings,
  key = "PC_",
  assay = "RNA"
)

message("\n=== Original Seurat Object ===")
print(seurat)
message("Cells: ", ncol(seurat))
message("Features: ", nrow(seurat))
message("Reductions: ", paste(names(seurat@reductions), collapse = ", "))

# Test output path
test_path <- file.path(getwd(), "test_seurat.scio")

message("\n=== Testing scio_write() with Seurat ===")
message("Writing to: ", test_path)
scio_write(seurat, path = test_path, overwrite = TRUE, compress = TRUE)

# Check what files were created
message("\n=== Checking created files ===")
if (dir.exists(test_path)) {
  files <- list.files(test_path, recursive = TRUE, full.names = FALSE)
  message("Files created:")
  for (f in files) {
    file_info <- file.info(file.path(test_path, f))
    message(sprintf("  %s (%.2f KB)", f, file_info$size / 1024))
  }
} else {
  message("ERROR: Output directory not created!")
}

message("\n=== Testing scio_read() back to Seurat ===")
message("Reading from: ", test_path)
seurat_back <- scio_read(test_path, output = "Seurat")

message("\n=== Loaded Seurat Object ===")
print(seurat_back)
message("Cells: ", ncol(seurat_back))
message("Features: ", nrow(seurat_back))
message("Reductions: ", paste(names(seurat_back@reductions), collapse = ", "))

# Verify data integrity
message("\n=== Verifying data integrity ===")
checks <- list(
  "Dimensions match" = all(dim(seurat) == dim(seurat_back)),
  "Count matrix matches" = all(GetAssayData(seurat, slot = "counts") == GetAssayData(seurat_back, slot = "counts")),
  "Metadata columns match" = all(colnames(seurat@meta.data) == colnames(seurat_back@meta.data)),
  "Cell names match" = all(colnames(seurat) == colnames(seurat_back)),
  "Gene names match" = all(rownames(seurat) == rownames(seurat_back)),
  "Reductions preserved" = "pca" %in% names(seurat_back@reductions)
)

for (check_name in names(checks)) {
  status <- if (checks[[check_name]]) "✓" else "✗"
  message(sprintf("%s %s", status, check_name))
}

# Clean up
message("\n=== Cleaning up ===")
unlink(test_path, recursive = TRUE)
message("Test output removed")

# Final summary
all_passed <- all(unlist(checks))
if (all_passed) {
  message("\n=== ✓ ALL SEURAT TESTS PASSED ===")
  quit(status = 0)
} else {
  message("\n=== ✗ SOME SEURAT TESTS FAILED ===")
  quit(status = 1)
}
