# Writers module - Save all components to folder structure
# Part of the scBridge R package

#' Save all components to folder structure
#'
#' Saves ALL 10 components:
#' X, obs, var, obsm, varm, obsp, varp, layers, raw, uns
#'
#' @param components Named list with all components
#' @param folder_path Output folder path
#' @param compress Whether to compress MTX files
#' @export
save_to_folder <- function(components, folder_path, compress = TRUE) {
  # Create output directory
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }

  saved_files <- list()

  # =========================================================================
  # 1. X - Expression matrix (transpose to genes × cells for storage)
  # =========================================================================
  matrix_file <- ifelse(compress, "matrix.mtx.gz", "matrix.mtx")
  X_transposed <- Matrix::t(components$X)  # Transpose to genes × cells
  save_mtx(X_transposed, file.path(folder_path, matrix_file), compress = compress)
  saved_files$X <- matrix_file

  # =========================================================================
  # 2. obs - Cell IDs (barcodes) + metadata
  # =========================================================================
  barcodes_file <- ifelse(compress, "barcodes.tsv.gz", "barcodes.tsv")
  barcodes_df <- data.frame(barcode = rownames(components$X))
  write.table(barcodes_df, file.path(folder_path, barcodes_file),
             sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  saved_files$barcodes <- barcodes_file

  # Save cell metadata as Parquet
  save_parquet(components$obs, file.path(folder_path, "obs.parquet"))
  saved_files$obs <- "obs.parquet"

  # =========================================================================
  # 3. var - Gene IDs (features) + metadata
  # =========================================================================
  features_file <- ifelse(compress, "features.tsv.gz", "features.tsv")
  features_df <- data.frame(
    gene_id = colnames(components$X),
    gene_name = colnames(components$X),
    feature_type = "Gene Expression"
  )
  write.table(features_df, file.path(folder_path, features_file),
             sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  saved_files$features <- features_file

  # Save gene metadata as Parquet
  if (ncol(components$var) > 0) {
    save_parquet(components$var, file.path(folder_path, "var.parquet"))
    saved_files$var <- "var.parquet"
  }

  # =========================================================================
  # 4. obsm - Cell embeddings
  # =========================================================================
  if (length(components$obsm) > 0) {
    obsm_dir <- file.path(folder_path, "obsm")
    if (!dir.exists(obsm_dir)) dir.create(obsm_dir)

    for (key in names(components$obsm)) {
      emb_df <- as.data.frame(components$obsm[[key]])
      rownames(emb_df) <- rownames(components$X)
      save_parquet(emb_df, file.path(obsm_dir, paste0(key, ".parquet")))
      saved_files[[paste0("obsm_", key)]] <- file.path("obsm", paste0(key, ".parquet"))
    }
  }

  # =========================================================================
  # 5. varm - Gene embeddings
  # =========================================================================
  if (length(components$varm) > 0) {
    varm_dir <- file.path(folder_path, "varm")
    if (!dir.exists(varm_dir)) dir.create(varm_dir)

    for (key in names(components$varm)) {
      varm_df <- as.data.frame(components$varm[[key]])
      rownames(varm_df) <- colnames(components$X)
      save_parquet(varm_df, file.path(varm_dir, paste0(key, ".parquet")))
      saved_files[[paste0("varm_", key)]] <- file.path("varm", paste0(key, ".parquet"))
    }
  }

  # =========================================================================
  # 6. obsp - Cell-cell graphs
  # =========================================================================
  if (length(components$obsp) > 0) {
    obsp_dir <- file.path(folder_path, "obsp")
    if (!dir.exists(obsp_dir)) dir.create(obsp_dir)

    for (key in names(components$obsp)) {
      graph_file <- ifelse(compress, paste0(key, ".mtx.gz"), paste0(key, ".mtx"))
      save_mtx(components$obsp[[key]], file.path(obsp_dir, graph_file), compress = compress)
      saved_files[[paste0("obsp_", key)]] <- file.path("obsp", graph_file)
    }
  }

  # =========================================================================
  # 7. varp - Gene-gene graphs
  # =========================================================================
  if (length(components$varp) > 0) {
    varp_dir <- file.path(folder_path, "varp")
    if (!dir.exists(varp_dir)) dir.create(varp_dir)

    for (key in names(components$varp)) {
      graph_file <- ifelse(compress, paste0(key, ".mtx.gz"), paste0(key, ".mtx"))
      save_mtx(components$varp[[key]], file.path(varp_dir, graph_file), compress = compress)
      saved_files[[paste0("varp_", key)]] <- file.path("varp", graph_file)
    }
  }

  # =========================================================================
  # 8. layers - Additional matrices
  # =========================================================================
  if (length(components$layers) > 0) {
    layers_dir <- file.path(folder_path, "layers")
    if (!dir.exists(layers_dir)) dir.create(layers_dir)

    for (key in names(components$layers)) {
      layer_file <- ifelse(compress, paste0(key, ".mtx.gz"), paste0(key, ".mtx"))
      layer_transposed <- Matrix::t(components$layers[[key]])  # Transpose to genes × cells
      save_mtx(layer_transposed, file.path(layers_dir, layer_file), compress = compress)
      saved_files[[paste0("layer_", key)]] <- file.path("layers", layer_file)
    }
  }

  # =========================================================================
  # 9. raw - Raw data
  # =========================================================================
  if (!is.null(components$raw)) {
    raw_dir <- file.path(folder_path, "raw")
    if (!dir.exists(raw_dir)) dir.create(raw_dir)

    # Save raw X matrix
    raw_matrix_file <- ifelse(compress, "matrix.mtx.gz", "matrix.mtx")
    raw_X_transposed <- Matrix::t(components$raw$X)  # Transpose to genes × cells
    save_mtx(raw_X_transposed, file.path(raw_dir, raw_matrix_file), compress = compress)
    saved_files$raw_X <- file.path("raw", raw_matrix_file)

    # Save raw var metadata
    if (ncol(components$raw$var) > 0) {
      save_parquet(components$raw$var, file.path(raw_dir, "var.parquet"))
      saved_files$raw_var <- file.path("raw", "var.parquet")
    }

    # Save raw gene IDs
    raw_features_file <- ifelse(compress, "features.tsv.gz", "features.tsv")
    raw_features_df <- data.frame(
      gene_id = colnames(components$raw$X),
      gene_name = colnames(components$raw$X),
      feature_type = "Gene Expression"
    )
    write.table(raw_features_df, file.path(raw_dir, raw_features_file),
               sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    saved_files$raw_features <- file.path("raw", raw_features_file)
  }

  # =========================================================================
  # 10. uns - Unstructured metadata
  # =========================================================================
  if (length(components$uns) > 0) {
    save_json(components$uns, file.path(folder_path, "uns.json"))
    saved_files$uns <- "uns.json"
  }

  # =========================================================================
  # Create manifest file
  # =========================================================================
  manifest <- list(
    format = "scBridge v1.0",
    created_by = "scBridge::write",
    dimensions = list(
      n_obs = nrow(components$X),
      n_vars = ncol(components$X)
    ),
    components = list(
      X = TRUE,
      obs = TRUE,
      var = ncol(components$var) > 0,
      obsm = names(components$obsm),
      varm = names(components$varm),
      obsp = names(components$obsp),
      varp = names(components$varp),
      layers = names(components$layers),
      raw = !is.null(components$raw),
      uns = names(components$uns)
    ),
    files = saved_files
  )

  save_json(manifest, file.path(folder_path, "manifest.json"))

  return(saved_files)
}
