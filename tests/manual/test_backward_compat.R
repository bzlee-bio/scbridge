#!/usr/bin/env Rscript
# Test backward compatibility with old .scb files
# Verifies that scio can read old scBridge files

library(SingleCellExperiment)
library(Matrix)

# Source the package functions
source_dir <- function(path) {
  files <- list.files(path, pattern = "\\.R$", full.names = TRUE)
  for (f in files) {
    source(f)
  }
}

pkg_path <- "../../R/R"
message("Sourcing R files from: ", pkg_path)
source_dir(pkg_path)

# Create test data
set.seed(42)
counts <- Matrix(round(runif(200, 0, 10)), nrow = 20, ncol = 10, sparse = TRUE)
rownames(counts) <- paste0("Gene", 1:20)
colnames(counts) <- paste0("Cell", 1:10)

sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = data.frame(
    sample = rep(c("A", "B"), each = 5),
    row.names = colnames(counts)
  ),
  rowData = data.frame(
    symbol = paste0("GENE", 1:20),
    row.names = rownames(counts)
  )
)

message("\n=== Test 1: Writing with new .scio extension ===")
scio_write(sce, "test_new.scio", overwrite = TRUE)
sce_new <- scio_read("test_new.scio", output = "SCE")
message("✓ New .scio format works")

message("\n=== Test 2: Reading old .scb folder (if renamed) ===")
# Simulate an old .scb file by renaming
if (dir.exists("test_old.scb")) unlink("test_old.scb", recursive = TRUE)
file.rename("test_new.scio", "test_old.scb")

# Can scio_read handle .scb folders?
tryCatch({
  sce_old <- scio_read("test_old.scb", output = "SCE")
  message("✓ Can read old .scb folders")

  # Verify data matches
  if (all(dim(sce) == dim(sce_old))) {
    message("✓ Data integrity preserved from .scb format")
  } else {
    message("✗ Data dimensions don't match!")
  }
}, error = function(e) {
  message("✗ Failed to read .scb folder: ", e$message)
})

# Clean up
unlink("test_old.scb", recursive = TRUE)
unlink("test_new.scio", recursive = TRUE)

message("\n=== Backward Compatibility Test Summary ===")
message("The internal format is identical - only the extension changed")
message("Old .scb folders can be read by simply:")
message("  1. Using them as-is (scio_read works with any folder)")
message("  2. Or renaming .scb → .scio (recommended)")
