# Readers module - Load all components from folder structure
# Part of the scio R package


#' Internal: Load sparse matrix with format auto-detection
#'
#' Detects v2.x binary CSC format vs v1.x MTX format and loads accordingly.
#' Handles orientation metadata from manifest for proper transpose behavior.
#'
#' @param folder_path Path to data folder
#' @param file_base Base file path from manifest (e.g., "matrix" or "matrix.mtx.gz")
#' @param manifest Manifest object
#' @param is_expression_matrix If TRUE, this is X/layer/raw.X that should be cells×genes
#' @param transpose For non-expression matrices (obsp/varp), whether to transpose
#' @return Sparse matrix
.load_sparse_matrix <- function(folder_path, file_base, manifest,
                                is_expression_matrix = FALSE, transpose = FALSE) {
  file_path <- file.path(folder_path, file_base)

  # Check orientation from manifest (v0.1.2+)
  orientation <- manifest$orientation
  if (is.null(orientation)) {
    orientation <- "genes_x_cells"  # Default for legacy formats
  }

  # Determine if we need to transpose
  if (is_expression_matrix) {
    # For expression matrices (X, layers, raw.X)
    if (orientation == "cells_x_genes") {
      # v0.1.2: already in cells×genes format - no transpose needed
      need_transpose <- FALSE
    } else {
      # v2.0 or v1.x: genes×cells format - need transpose
      need_transpose <- TRUE
    }
  } else {
    # For non-expression matrices (obsp, varp), use the provided transpose flag
    need_transpose <- transpose
  }

  # Try v2.x binary CSC format first
  shape_file <- paste0(file_path, ".shape.json")
  if (file.exists(shape_file)) {
    # v2.x binary CSC format
    return(load_sparse_binary(file_path, transpose = need_transpose))
  }

  # Fallback to v1.x MTX format
  mtx_file <- file_path
  if (!file.exists(mtx_file)) {
    # Try with .gz extension
    mtx_file <- paste0(file_path, ".gz")
  }

  mat <- load_mtx(mtx_file)
  if (need_transpose) {
    mat <- Matrix::t(mat)
  }
  return(mat)
}


#' Load all data components from folder structure
#'
#' This loads ALL 10 components:
#' X, obs, var, obsm, varm, obsp, varp, layers, raw, uns
#'
#' @param folder_path Path to folder containing saved data
#' @return Named list with all components
#' @export
load_from_folder <- function(folder_path) {
  if (!dir.exists(folder_path)) {
    stop(paste("Folder not found:", folder_path))
  }

  # Load manifest
  manifest_file <- file.path(folder_path, "manifest.json")
  if (!file.exists(manifest_file)) {
    stop(paste("Manifest file not found:", manifest_file))
  }

  manifest <- load_json(manifest_file)

  result <- list()

  # =========================================================================
  # 1. Load X matrix (cells × genes)
  # v0.1.2: already cells×genes, no transpose needed
  # v2.0/v1.x: genes×cells, transpose needed (handled by is_expression_matrix flag)
  # =========================================================================
  result$X <- .load_sparse_matrix(folder_path, manifest$files$X, manifest,
                                   is_expression_matrix = TRUE)

  # =========================================================================
  # 2. Load cell IDs (obs_names)
  # =========================================================================
  barcodes_file <- file.path(folder_path, manifest$files$barcodes)
  barcodes_df <- read.table(barcodes_file, header = FALSE, stringsAsFactors = FALSE)
  obs_names <- barcodes_df$V1
  result$obs_names <- obs_names

  # =========================================================================
  # 3. Load gene IDs (var_names)
  # =========================================================================
  features_file <- file.path(folder_path, manifest$files$features)
  features_df <- read.table(features_file, header = FALSE, sep = "\t",
                           stringsAsFactors = FALSE)
  var_names <- features_df$V1
  result$var_names <- var_names

  # Set row and column names for X
  rownames(result$X) <- obs_names
  colnames(result$X) <- var_names

  # =========================================================================
  # 4. Load obs (cell metadata)
  # =========================================================================
  obs_file <- file.path(folder_path, manifest$files$obs)
  result$obs <- load_parquet(obs_file)
  rownames(result$obs) <- obs_names

  # =========================================================================
  # 5. Load var (gene metadata)
  # =========================================================================
  if (!is.null(manifest$files$var)) {
    var_file <- file.path(folder_path, manifest$files$var)
    result$var <- load_parquet(var_file)
    rownames(result$var) <- var_names
  } else {
    result$var <- data.frame(row.names = var_names)
  }

  # =========================================================================
  # 6. Load obsm (cell embeddings)
  # =========================================================================
  result$obsm <- list()
  if (length(manifest$components$obsm) > 0) {
    for (key in manifest$components$obsm) {
      file_key <- paste0("obsm_", key)
      file_path <- file.path(folder_path, manifest$files[[file_key]])
      emb_df <- load_parquet(file_path)
      result$obsm[[key]] <- as.matrix(emb_df)
      rownames(result$obsm[[key]]) <- obs_names
    }
  }

  # =========================================================================
  # 7. Load varm (gene embeddings)
  # =========================================================================
  result$varm <- list()
  if (length(manifest$components$varm) > 0) {
    for (key in manifest$components$varm) {
      file_key <- paste0("varm_", key)
      file_path <- file.path(folder_path, manifest$files[[file_key]])
      varm_df <- load_parquet(file_path)
      result$varm[[key]] <- as.matrix(varm_df)
      rownames(result$varm[[key]]) <- var_names
    }
  }

  # =========================================================================
  # 8. Load obsp (cell-cell graphs) - symmetric, no transpose needed
  # =========================================================================
  result$obsp <- list()
  if (length(manifest$components$obsp) > 0) {
    for (key in manifest$components$obsp) {
      file_key <- paste0("obsp_", key)
      result$obsp[[key]] <- .load_sparse_matrix(folder_path, manifest$files[[file_key]],
                                                 manifest, transpose = FALSE)
      rownames(result$obsp[[key]]) <- obs_names
      colnames(result$obsp[[key]]) <- obs_names
    }
  }

  # =========================================================================
  # 9. Load varp (gene-gene graphs) - symmetric, no transpose needed
  # =========================================================================
  result$varp <- list()
  if (length(manifest$components$varp) > 0) {
    for (key in manifest$components$varp) {
      file_key <- paste0("varp_", key)
      result$varp[[key]] <- .load_sparse_matrix(folder_path, manifest$files[[file_key]],
                                                 manifest, transpose = FALSE)
      rownames(result$varp[[key]]) <- var_names
      colnames(result$varp[[key]]) <- var_names
    }
  }

  # =========================================================================
  # 10. Load layers (additional matrices)
  # v0.1.2: already cells×genes, no transpose needed
  # v2.0/v1.x: genes×cells, transpose needed (handled by is_expression_matrix flag)
  # =========================================================================
  result$layers <- list()
  if (length(manifest$components$layers) > 0) {
    for (key in manifest$components$layers) {
      file_key <- paste0("layer_", key)
      result$layers[[key]] <- .load_sparse_matrix(folder_path, manifest$files[[file_key]],
                                                   manifest, is_expression_matrix = TRUE)
      rownames(result$layers[[key]]) <- obs_names
      colnames(result$layers[[key]]) <- var_names
    }
  }

  # =========================================================================
  # 11. Load raw data (if present)
  # v0.1.2: raw.X already cells×genes, no transpose needed
  # v2.0/v1.x: genes×cells, transpose needed (handled by is_expression_matrix flag)
  # =========================================================================
  if (manifest$components$raw) {
    result$raw <- list()

    # Load raw X matrix
    result$raw$X <- .load_sparse_matrix(folder_path, manifest$files$raw_X,
                                         manifest, is_expression_matrix = TRUE)

    # Load raw gene IDs
    raw_features_file <- file.path(folder_path, strsplit(manifest$files$raw_features, "/")[[1]][2])
    raw_features_path <- file.path(folder_path, "raw", basename(raw_features_file))
    raw_features_df <- read.table(raw_features_path, header = FALSE, sep = "\t",
                                  stringsAsFactors = FALSE)
    raw_var_names <- raw_features_df$V1

    # Load raw var metadata if present
    if (!is.null(manifest$files$raw_var)) {
      raw_var_file <- file.path(folder_path, manifest$files$raw_var)
      result$raw$var <- load_parquet(raw_var_file)
      rownames(result$raw$var) <- raw_var_names
    } else {
      result$raw$var <- data.frame(row.names = raw_var_names)
    }

    # Set row and column names for raw X
    rownames(result$raw$X) <- obs_names
    colnames(result$raw$X) <- raw_var_names
  }

  # =========================================================================
  # 12. Load uns (unstructured metadata)
  # =========================================================================
  if (!is.null(manifest$files$uns)) {
    uns_file <- file.path(folder_path, manifest$files$uns)
    result$uns <- load_json(uns_file)
  } else {
    result$uns <- list()
  }

  return(result)
}
