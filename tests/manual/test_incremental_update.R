#!/usr/bin/env Rscript
# Test incremental update feature for scio
# Tests hash-based change detection and selective component updates

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
pkg_path <- "../../R/R"
message("Sourcing R files from: ", pkg_path)
source_dir(pkg_path)

# Create a test SCE object
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
    row.names = rownames(counts)
  )
)

# Add PCA embedding
pca <- matrix(rnorm(10 * 5), ncol = 5)
rownames(pca) <- colnames(sce)
colnames(pca) <- paste0("PC", 1:5)
reducedDim(sce, "PCA") <- pca

test_path <- file.path(getwd(), "test_incremental.scio")

# Clean up any existing test file
if (dir.exists(test_path)) unlink(test_path, recursive = TRUE)

message("\n=== Test 1: Initial save with hash computation ===")
scio_write(sce, test_path)

# Check manifest has hashes
manifest <- jsonlite::fromJSON(file.path(test_path, "manifest.json"))
if (!is.null(manifest$hashes)) {
  message("  Hashes found in manifest:")
  for (key in names(manifest$hashes)) {
    message("    ", key, ": ", substr(manifest$hashes[[key]], 1, 8), "...")
  }
} else {
  message("  ERROR: No hashes in manifest!")
}

message("\n=== Test 2: Update with no changes (should skip all) ===")
scio_write(sce, test_path, update = TRUE)

message("\n=== Test 3: Update obs only ===")
# Modify cell metadata
sce$new_column <- 1:10
scio_write(sce, test_path, update = TRUE)

# Verify only obs was updated
manifest2 <- jsonlite::fromJSON(file.path(test_path, "manifest.json"))
message("  Checking obs hash changed: ", manifest$hashes$obs != manifest2$hashes$obs)
message("  Checking X hash unchanged: ", manifest$hashes$X == manifest2$hashes$X)

message("\n=== Test 4: Add new obsm (UMAP) ===")
umap <- matrix(rnorm(10 * 2), ncol = 2)
rownames(umap) <- colnames(sce)
colnames(umap) <- c("UMAP1", "UMAP2")
reducedDim(sce, "UMAP") <- umap

scio_write(sce, test_path, update = TRUE)

# Verify UMAP was added
manifest3 <- jsonlite::fromJSON(file.path(test_path, "manifest.json"))
message("  UMAP hash present: ", "obsm_UMAP" %in% names(manifest3$hashes))

message("\n=== Test 5: Verify data integrity after updates ===")
sce_back <- scio_read(test_path, output = "SCE")

checks <- list(
  "Dimensions match" = all(dim(sce) == dim(sce_back)),
  "new_column preserved" = "new_column" %in% colnames(colData(sce_back)),
  "PCA preserved" = "PCA" %in% reducedDimNames(sce_back),
  "UMAP preserved" = "UMAP" %in% reducedDimNames(sce_back),
  "PCA values match" = all(abs(reducedDim(sce, "PCA") - reducedDim(sce_back, "PCA")) < 1e-10),
  "UMAP values match" = all(abs(reducedDim(sce, "UMAP") - reducedDim(sce_back, "UMAP")) < 1e-10)
)

for (check_name in names(checks)) {
  status <- if (checks[[check_name]]) "OK" else "FAIL"
  message(sprintf("  %s %s", status, check_name))
}

message("\n=== Test 6: Test overwrite=TRUE vs update=TRUE conflict ===")
tryCatch({
  scio_write(sce, test_path, overwrite = TRUE, update = TRUE)
  message("  ERROR: Should have thrown an error!")
}, error = function(e) {
  message("  Correctly caught error: ", e$message)
})

message("\n=== Test 7: Test update on non-existent file ===")
nonexistent_path <- file.path(getwd(), "nonexistent.scio")
tryCatch({
  scio_write(sce, nonexistent_path, update = TRUE)
  message("  New file created with update=TRUE (first save)")
  if (dir.exists(nonexistent_path)) {
    unlink(nonexistent_path, recursive = TRUE)
  }
}, error = function(e) {
  message("  ERROR: ", e$message)
})

# Clean up
message("\n=== Cleaning up ===")
unlink(test_path, recursive = TRUE)
message("Test output removed")

# Final summary
all_passed <- all(unlist(checks))
if (all_passed) {
  message("\n=== ALL INCREMENTAL UPDATE TESTS PASSED ===")
  quit(status = 0)
} else {
  message("\n=== SOME INCREMENTAL UPDATE TESTS FAILED ===")
  quit(status = 1)
}
