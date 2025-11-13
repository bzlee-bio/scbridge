# Readers module - Load all components from folder structure
# Part of the scBridge R package

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
  # 1. Load X matrix (transpose from genes × cells to cells × genes)
  # =========================================================================
  matrix_file <- file.path(folder_path, manifest$files$X)
  X <- load_mtx(matrix_file)
  X <- Matrix::t(X)  # Transpose to cells × genes
  result$X <- X

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
  # 8. Load obsp (cell-cell graphs)
  # =========================================================================
  result$obsp <- list()
  if (length(manifest$components$obsp) > 0) {
    for (key in manifest$components$obsp) {
      file_key <- paste0("obsp_", key)
      file_path <- file.path(folder_path, manifest$files[[file_key]])
      result$obsp[[key]] <- load_mtx(file_path)
      rownames(result$obsp[[key]]) <- obs_names
      colnames(result$obsp[[key]]) <- obs_names
    }
  }

  # =========================================================================
  # 9. Load varp (gene-gene graphs)
  # =========================================================================
  result$varp <- list()
  if (length(manifest$components$varp) > 0) {
    for (key in manifest$components$varp) {
      file_key <- paste0("varp_", key)
      file_path <- file.path(folder_path, manifest$files[[file_key]])
      result$varp[[key]] <- load_mtx(file_path)
      rownames(result$varp[[key]]) <- var_names
      colnames(result$varp[[key]]) <- var_names
    }
  }

  # =========================================================================
  # 10. Load layers (additional matrices)
  # =========================================================================
  result$layers <- list()
  if (length(manifest$components$layers) > 0) {
    for (key in manifest$components$layers) {
      file_key <- paste0("layer_", key)
      file_path <- file.path(folder_path, manifest$files[[file_key]])
      layer_mat <- load_mtx(file_path)
      result$layers[[key]] <- Matrix::t(layer_mat)  # Transpose to cells × genes
      rownames(result$layers[[key]]) <- obs_names
      colnames(result$layers[[key]]) <- var_names
    }
  }

  # =========================================================================
  # 11. Load raw data (if present)
  # =========================================================================
  if (manifest$components$raw) {
    result$raw <- list()

    # Load raw X matrix
    raw_matrix_file <- file.path(folder_path, strsplit(manifest$files$raw_X, "/")[[1]][2])
    raw_matrix_path <- file.path(folder_path, "raw", basename(raw_matrix_file))
    raw_X <- load_mtx(raw_matrix_path)
    result$raw$X <- Matrix::t(raw_X)  # Transpose to cells × genes

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
