#!/usr/bin/env Rscript
# Manual test for scio write/read round-trip
# Tests that scio_write() and scio_read() work correctly with compression

library(SingleCellExperiment)
library(Matrix)

# Source the package functions directly (for development testing)
source_dir <- function(path) {
  files <- list.files(path, pattern = "\\.R$", full.names = TRUE)
  for (f in files) {
    source(f)
  }
}

# Source all R files from the package
pkg_path <- file.path(dirname(dirname(dirname(getwd()))), "R", "R")
if (!dir.exists(pkg_path)) {
  pkg_path <- "../../R/R"  # Try relative path
}
message("Sourcing R files from: ", pkg_path)
source_dir(pkg_path)

# Create a small test SCE object
set.seed(42)
counts <- Matrix(round(runif(200, 0, 10)), nrow = 20, ncol = 10, sparse = TRUE)
rownames(counts) <- paste0("Gene", 1:20)
colnames(counts) <- paste0("Cell", 1:10)

sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = data.frame(
    sample = rep(c("A", "B"), each = 5),
    celltype = rep(c("T", "B", "NK"), length.out = 10),
    row.names = colnames(counts)
  ),
  rowData = data.frame(
    symbol = paste0("GENE", 1:20),
    type = "Gene Expression",
    row.names = rownames(counts)
  )
)

# Add some reduced dimensions
reducedDim(sce, "PCA") <- matrix(rnorm(10 * 5), ncol = 5)
colnames(reducedDim(sce, "PCA")) <- paste0("PC", 1:5)

message("\n=== Original SCE Object ===")
print(sce)
message("Dims: ", nrow(sce), " x ", ncol(sce))
message("Assays: ", paste(assayNames(sce), collapse = ", "))
message("ReducedDims: ", paste(reducedDimNames(sce), collapse = ", "))

# Test output path
test_path <- file.path(getwd(), "test_output.scio")

message("\n=== Testing scio_write() ===")
message("Writing to: ", test_path)
scio_write(sce, path = test_path, overwrite = TRUE, compress = TRUE)

# Check what files were created
message("\n=== Checking created files ===")
if (dir.exists(test_path)) {
  files <- list.files(test_path, recursive = TRUE, full.names = FALSE)
  message("Files created:")
  for (f in files) {
    file_info <- file.info(file.path(test_path, f))
    message(sprintf("  %s (%.2f KB)", f, file_info$size / 1024))
  }

  # Check for any numbered files (the bug we're fixing)
  numbered_files <- grep("^[0-9]+$", files, value = TRUE)
  if (length(numbered_files) > 0) {
    message("\n!!! ERROR: Found numbered files (bug not fixed): ",
            paste(numbered_files, collapse = ", "))
  } else {
    message("\n✓ No numbered files found (bug is fixed!)")
  }
} else {
  message("ERROR: Output directory not created!")
}

message("\n=== Testing scio_read() ===")
message("Reading from: ", test_path)
sce_back <- scio_read(test_path, output = "SCE")

message("\n=== Loaded SCE Object ===")
print(sce_back)
message("Dims: ", nrow(sce_back), " x ", ncol(sce_back))
message("Assays: ", paste(assayNames(sce_back), collapse = ", "))
message("ReducedDims: ", paste(reducedDimNames(sce_back), collapse = ", "))

# Verify data integrity
message("\n=== Verifying data integrity ===")
checks <- list(
  "Dimensions match" = all(dim(sce) == dim(sce_back)),
  "Count matrix matches" = all(assay(sce, "counts") == assay(sce_back, "counts")),
  "Cell metadata matches" = all(colData(sce) == colData(sce_back)),
  "Gene metadata matches" = nrow(rowData(sce)) == nrow(rowData(sce_back)),
  "PCA embeddings match" = all(abs(reducedDim(sce, "PCA") - reducedDim(sce_back, "PCA")) < 1e-10)
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
  message("\n=== ✓ ALL TESTS PASSED ===")
  quit(status = 0)
} else {
  message("\n=== ✗ SOME TESTS FAILED ===")
  quit(status = 1)
}
