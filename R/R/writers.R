# Writers module - Save all components to folder structure
# Part of the scio R package

#' Save all components to folder structure
#'
#' Saves ALL 10 components using binary sparse format (v0.1.3):
#' X, obs, var, obsm, varm, obsp, varp, layers, raw, uns
#'
#' @param components Named list with all components
#' @param folder_path Output folder path
#' @param sparse_format Format for sparse matrices: "csr" (default) or "csc"
#' @param compute_hashes Whether to compute and store MD5 hashes for
#'   incremental updates (default: FALSE for backward compatibility)
#' @export
save_to_folder <- function(components, folder_path, sparse_format = "csr", compute_hashes = FALSE) {
  # Create output directory
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }

  saved_files <- list()
  hashes <- list()  # Store hashes for incremental updates

  # =========================================================================
  # 1. X - Expression matrix (cells × genes orientation, binary format)
  # =========================================================================
  # Save in cells × genes orientation (no transpose needed for Python read)
  save_sparse_binary(components$X, file.path(folder_path, "matrix"),
                     sparse_format = sparse_format, orientation = "cells_x_genes")
  saved_files$X <- "matrix"
  if (compute_hashes) hashes$X <- compute_hash(components$X)

  # =========================================================================
  # 2. obs - Cell IDs (barcodes) + metadata
  # =========================================================================
  # Always gzip barcodes - small file, fast compression
  barcodes_file <- "barcodes.tsv.gz"
  barcodes_df <- data.frame(barcode = rownames(components$X))
  con <- gzfile(file.path(folder_path, barcodes_file), "wb")
  write.table(barcodes_df, con, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  close(con)
  saved_files$barcodes <- barcodes_file

  # Save cell metadata as Parquet
  save_parquet(components$obs, file.path(folder_path, "obs.parquet"))
  saved_files$obs <- "obs.parquet"
  if (compute_hashes) hashes$obs <- compute_hash(components$obs)

  # =========================================================================
  # 3. var - Gene IDs (features) + metadata
  # =========================================================================
  # Always gzip features - small file, fast compression
  features_file <- "features.tsv.gz"
  features_df <- data.frame(
    gene_id = colnames(components$X),
    gene_name = colnames(components$X),
    feature_type = "Gene Expression"
  )
  con <- gzfile(file.path(folder_path, features_file), "wb")
  write.table(features_df, con, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  close(con)
  saved_files$features <- features_file

  # Save gene metadata as Parquet
  if (ncol(components$var) > 0) {
    save_parquet(components$var, file.path(folder_path, "var.parquet"))
    saved_files$var <- "var.parquet"
    if (compute_hashes) hashes$var <- compute_hash(components$var)
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
      if (compute_hashes) hashes[[paste0("obsm_", key)]] <- compute_hash(components$obsm[[key]])
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
      if (compute_hashes) hashes[[paste0("varm_", key)]] <- compute_hash(components$varm[[key]])
    }
  }

  # =========================================================================
  # 6. obsp - Cell-cell graphs (binary format)
  # =========================================================================
  if (length(components$obsp) > 0) {
    obsp_dir <- file.path(folder_path, "obsp")
    if (!dir.exists(obsp_dir)) dir.create(obsp_dir)

    for (key in names(components$obsp)) {
      save_sparse_binary(components$obsp[[key]], file.path(obsp_dir, key),
                         sparse_format = sparse_format)
      saved_files[[paste0("obsp_", key)]] <- file.path("obsp", key)
      if (compute_hashes) hashes[[paste0("obsp_", key)]] <- compute_hash(components$obsp[[key]])
    }
  }

  # =========================================================================
  # 7. varp - Gene-gene graphs (binary format)
  # =========================================================================
  if (length(components$varp) > 0) {
    varp_dir <- file.path(folder_path, "varp")
    if (!dir.exists(varp_dir)) dir.create(varp_dir)

    for (key in names(components$varp)) {
      save_sparse_binary(components$varp[[key]], file.path(varp_dir, key),
                         sparse_format = sparse_format)
      saved_files[[paste0("varp_", key)]] <- file.path("varp", key)
      if (compute_hashes) hashes[[paste0("varp_", key)]] <- compute_hash(components$varp[[key]])
    }
  }

  # =========================================================================
  # 8. layers - Additional matrices (binary format, cells × genes)
  # =========================================================================
  if (length(components$layers) > 0) {
    layers_dir <- file.path(folder_path, "layers")
    if (!dir.exists(layers_dir)) dir.create(layers_dir)

    for (key in names(components$layers)) {
      # Save in cells × genes orientation (same as X matrix)
      save_sparse_binary(components$layers[[key]], file.path(layers_dir, key),
                         sparse_format = sparse_format, orientation = "cells_x_genes")
      saved_files[[paste0("layer_", key)]] <- file.path("layers", key)
      if (compute_hashes) hashes[[paste0("layer_", key)]] <- compute_hash(components$layers[[key]])
    }
  }

  # =========================================================================
  # 9. raw - Raw data (binary format, cells × genes)
  # =========================================================================
  if (!is.null(components$raw)) {
    raw_dir <- file.path(folder_path, "raw")
    if (!dir.exists(raw_dir)) dir.create(raw_dir)

    # Save raw X matrix in cells × genes orientation
    save_sparse_binary(components$raw$X, file.path(raw_dir, "matrix"),
                       sparse_format = sparse_format, orientation = "cells_x_genes")
    saved_files$raw_X <- file.path("raw", "matrix")
    if (compute_hashes) hashes$raw_X <- compute_hash(components$raw$X)

    # Save raw var metadata
    if (ncol(components$raw$var) > 0) {
      save_parquet(components$raw$var, file.path(raw_dir, "var.parquet"))
      saved_files$raw_var <- file.path("raw", "var.parquet")
      if (compute_hashes) hashes$raw_var <- compute_hash(components$raw$var)
    }

    # Save raw gene IDs (always gzip - small file, fast compression)
    raw_features_file <- "features.tsv.gz"
    raw_features_df <- data.frame(
      gene_id = colnames(components$raw$X),
      gene_name = colnames(components$raw$X),
      feature_type = "Gene Expression"
    )
    con <- gzfile(file.path(raw_dir, raw_features_file), "wb")
    write.table(raw_features_df, con, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    close(con)
    saved_files$raw_features <- file.path("raw", raw_features_file)
  }

  # =========================================================================
  # 10. uns - Unstructured metadata
  # =========================================================================
  if (length(components$uns) > 0) {
    save_json(components$uns, file.path(folder_path, "uns.json"))
    saved_files$uns <- "uns.json"
    if (compute_hashes) hashes$uns <- compute_hash(components$uns)
  }

  # =========================================================================
  # Create manifest file
  # =========================================================================
  # Helper to ensure arrays are preserved as arrays in JSON (not unboxed to scalars)
  as_json_array <- function(x) {
    if (length(x) == 0) return(list())
    I(x)  # I() prevents auto_unbox from converting single elements to scalars
  }

  manifest <- list(
    format = "scio v0.1.3",
    created_by = "scio::write",
    orientation = "cells_x_genes",
    dimensions = list(
      n_obs = nrow(components$X),
      n_vars = ncol(components$X)
    ),
    components = list(
      X = TRUE,
      obs = TRUE,
      var = ncol(components$var) > 0,
      obsm = as_json_array(names(components$obsm)),
      varm = as_json_array(names(components$varm)),
      obsp = as_json_array(names(components$obsp)),
      varp = as_json_array(names(components$varp)),
      layers = as_json_array(names(components$layers)),
      raw = !is.null(components$raw),
      uns = as_json_array(names(components$uns))
    ),
    files = saved_files
  )

  # Add hashes if computed
  if (compute_hashes && length(hashes) > 0) {
    manifest$hashes <- hashes
  }

  save_json(manifest, file.path(folder_path, "manifest.json"))

  return(saved_files)
}


#' Update folder with only changed components (incremental update)
#'
#' Compares current components with stored hashes and only rewrites
#' components that have changed. Uses binary sparse format (v0.1.3).
#'
#' @param components Named list with all components
#' @param folder_path Existing .scio folder path
#' @param manifest Existing manifest with hashes
#' @param sparse_format Format for sparse matrices: "csr" (default) or "csc"
#' @export
update_folder <- function(components, folder_path, manifest, sparse_format = "csr") {
  old_hashes <- manifest$hashes
  new_hashes <- list()
  saved_files <- manifest$files  # Keep existing file paths
  updated_count <- 0
  skipped_count <- 0

  # Helper function to check if component changed
  check_and_update <- function(component, key, save_func) {
    new_hash <- compute_hash(component)
    new_hashes[[key]] <<- new_hash

    if (is.null(old_hashes[[key]]) || old_hashes[[key]] != new_hash) {
      save_func()
      updated_count <<- updated_count + 1
      message("    Updated: ", key)
      return(TRUE)
    } else {
      skipped_count <<- skipped_count + 1
      return(FALSE)
    }
  }

  # =========================================================================
  # 1. X - Expression matrix (binary format, cells × genes)
  # =========================================================================
  check_and_update(components$X, "X", function() {
    save_sparse_binary(components$X, file.path(folder_path, "matrix"),
                       sparse_format = sparse_format, orientation = "cells_x_genes")
  })

  # =========================================================================
  # 2. obs - Cell metadata
  # =========================================================================
  check_and_update(components$obs, "obs", function() {
    save_parquet(components$obs, file.path(folder_path, "obs.parquet"))
    # Also update barcodes (always gzip - small file)
    barcodes_file <- "barcodes.tsv.gz"
    barcodes_df <- data.frame(barcode = rownames(components$X))
    con <- gzfile(file.path(folder_path, barcodes_file), "wb")
    write.table(barcodes_df, con, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    close(con)
  })

  # =========================================================================
  # 3. var - Gene metadata
  # =========================================================================
  if (ncol(components$var) > 0) {
    check_and_update(components$var, "var", function() {
      save_parquet(components$var, file.path(folder_path, "var.parquet"))
    })
  }

  # =========================================================================
  # 4. obsm - Cell embeddings
  # =========================================================================
  if (length(components$obsm) > 0) {
    obsm_dir <- file.path(folder_path, "obsm")
    if (!dir.exists(obsm_dir)) dir.create(obsm_dir)

    for (key in names(components$obsm)) {
      hash_key <- paste0("obsm_", key)
      file_path <- file.path("obsm", paste0(key, ".parquet"))
      saved_files[[hash_key]] <- file_path  # Track file path
      check_and_update(components$obsm[[key]], hash_key, function() {
        emb_df <- as.data.frame(components$obsm[[key]])
        rownames(emb_df) <- rownames(components$X)
        save_parquet(emb_df, file.path(obsm_dir, paste0(key, ".parquet")))
      })
    }
  }

  # =========================================================================
  # 5. varm - Gene embeddings
  # =========================================================================
  if (length(components$varm) > 0) {
    varm_dir <- file.path(folder_path, "varm")
    if (!dir.exists(varm_dir)) dir.create(varm_dir)

    for (key in names(components$varm)) {
      hash_key <- paste0("varm_", key)
      file_path <- file.path("varm", paste0(key, ".parquet"))
      saved_files[[hash_key]] <- file_path
      check_and_update(components$varm[[key]], hash_key, function() {
        varm_df <- as.data.frame(components$varm[[key]])
        rownames(varm_df) <- colnames(components$X)
        save_parquet(varm_df, file.path(varm_dir, paste0(key, ".parquet")))
      })
    }
  }

  # =========================================================================
  # 6. obsp - Cell-cell graphs (binary format)
  # =========================================================================
  if (length(components$obsp) > 0) {
    obsp_dir <- file.path(folder_path, "obsp")
    if (!dir.exists(obsp_dir)) dir.create(obsp_dir)

    for (key in names(components$obsp)) {
      hash_key <- paste0("obsp_", key)
      saved_files[[hash_key]] <- file.path("obsp", key)
      check_and_update(components$obsp[[key]], hash_key, function() {
        save_sparse_binary(components$obsp[[key]], file.path(obsp_dir, key),
                           sparse_format = sparse_format)
      })
    }
  }

  # =========================================================================
  # 7. varp - Gene-gene graphs (binary format)
  # =========================================================================
  if (length(components$varp) > 0) {
    varp_dir <- file.path(folder_path, "varp")
    if (!dir.exists(varp_dir)) dir.create(varp_dir)

    for (key in names(components$varp)) {
      hash_key <- paste0("varp_", key)
      saved_files[[hash_key]] <- file.path("varp", key)
      check_and_update(components$varp[[key]], hash_key, function() {
        save_sparse_binary(components$varp[[key]], file.path(varp_dir, key),
                           sparse_format = sparse_format)
      })
    }
  }

  # =========================================================================
  # 8. layers - Additional matrices (binary format, cells × genes)
  # =========================================================================
  if (length(components$layers) > 0) {
    layers_dir <- file.path(folder_path, "layers")
    if (!dir.exists(layers_dir)) dir.create(layers_dir)

    for (key in names(components$layers)) {
      hash_key <- paste0("layer_", key)
      saved_files[[hash_key]] <- file.path("layers", key)
      check_and_update(components$layers[[key]], hash_key, function() {
        save_sparse_binary(components$layers[[key]], file.path(layers_dir, key),
                           sparse_format = sparse_format, orientation = "cells_x_genes")
      })
    }
  }

  # =========================================================================
  # 9. raw - Raw data (binary format, cells × genes)
  # =========================================================================
  if (!is.null(components$raw)) {
    raw_dir <- file.path(folder_path, "raw")
    if (!dir.exists(raw_dir)) dir.create(raw_dir)

    check_and_update(components$raw$X, "raw_X", function() {
      save_sparse_binary(components$raw$X, file.path(raw_dir, "matrix"),
                         sparse_format = sparse_format, orientation = "cells_x_genes")
    })

    if (ncol(components$raw$var) > 0) {
      check_and_update(components$raw$var, "raw_var", function() {
        save_parquet(components$raw$var, file.path(raw_dir, "var.parquet"))
      })
    }
  }

  # =========================================================================
  # 10. uns - Unstructured metadata
  # =========================================================================
  if (length(components$uns) > 0) {
    check_and_update(components$uns, "uns", function() {
      save_json(components$uns, file.path(folder_path, "uns.json"))
    })
  }

  # =========================================================================
  # Update manifest with new hashes, files, and component info
  # =========================================================================
  # Helper to ensure arrays are preserved as arrays in JSON (not unboxed to scalars)
  as_json_array <- function(x) {
    if (length(x) == 0) return(list())
    I(x)  # I() prevents auto_unbox from converting single elements to scalars
  }

  manifest$format <- "scio v0.1.3"
  manifest$orientation <- "cells_x_genes"
  manifest$hashes <- new_hashes
  manifest$files <- saved_files  # Update file paths
  manifest$last_updated <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")

  # Update components list with current component names
  manifest$components$obsm <- as_json_array(names(components$obsm))
  manifest$components$varm <- as_json_array(names(components$varm))
  manifest$components$obsp <- as_json_array(names(components$obsp))
  manifest$components$varp <- as_json_array(names(components$varp))
  manifest$components$layers <- as_json_array(names(components$layers))
  manifest$components$raw <- !is.null(components$raw)
  manifest$components$uns <- as_json_array(names(components$uns))

  save_json(manifest, file.path(folder_path, "manifest.json"))

  message("    Summary: ", updated_count, " components updated, ", skipped_count, " unchanged")

  return(invisible(list(updated = updated_count, skipped = skipped_count)))
}
