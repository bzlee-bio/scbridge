# Seurat object converters
# Part of the scBridge R package

#' Convert components list to Seurat object
#'
#' @param components Named list with X, obs, var, obsm, obsp, etc.
#' @return Seurat object
#' @export
create_seurat_from_components <- function(components) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package required. Install with: install.packages('Seurat')")
  }

  # Create Seurat object from count matrix (transpose to genes × cells)
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = Matrix::t(components$X),
    meta.data = components$obs
  )

  # Add embeddings as DimReduc objects
  if (length(components$obsm) > 0) {
    for (emb_name in names(components$obsm)) {
      # Determine key for Seurat (e.g., X_pca -> PC_, X_umap -> UMAP_)
      if (startsWith(emb_name, "X_")) {
        key <- toupper(sub("^X_", "", emb_name))
        key <- paste0(key, "_")
      } else {
        key <- paste0(toupper(emb_name), "_")
      }

      seurat_obj[[emb_name]] <- Seurat::CreateDimReducObject(
        embeddings = components$obsm[[emb_name]],
        key = key
      )
    }
  }

  # Add graphs (obsp) to Seurat object
  if (length(components$obsp) > 0) {
    for (graph_name in names(components$obsp)) {
      seurat_obj[[graph_name]] <- components$obsp[[graph_name]]
    }
  }

  # Add layers if present
  if (length(components$layers) > 0) {
    for (layer_name in names(components$layers)) {
      # Add as assay layer (transpose to genes × cells)
      seurat_obj[[layer_name]] <- Seurat::CreateAssayObject(
        counts = Matrix::t(components$layers[[layer_name]])
      )
    }
  }

  # Add raw data if present
  if (!is.null(components$raw)) {
    # Store in misc slot for now
    seurat_obj@misc$raw <- components$raw
  }

  # Add uns metadata
  if (length(components$uns) > 0) {
    seurat_obj@misc$uns <- components$uns
  }

  return(seurat_obj)
}


#' Extract components from Seurat object
#'
#' @param seurat_obj Seurat object
#' @return Named list with all components
#' @export
extract_components_from_seurat <- function(seurat_obj) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package required.")
  }

  components <- list()

  # Get count matrix (transpose from genes × cells to cells × genes)
  default_assay <- Seurat::DefaultAssay(seurat_obj)
  components$X <- Matrix::t(Seurat::GetAssayData(seurat_obj, slot = "counts", assay = default_assay))

  # Get cell metadata
  components$obs <- seurat_obj@meta.data

  # Get gene metadata from features
  components$var <- data.frame(row.names = rownames(seurat_obj))

  # Get embeddings (obsm)
  components$obsm <- list()
  reductions <- names(seurat_obj@reductions)
  if (length(reductions) > 0) {
    for (reduction in reductions) {
      components$obsm[[reduction]] <- Seurat::Embeddings(seurat_obj, reduction = reduction)
    }
  }

  # Get graphs (obsp)
  components$obsp <- list()
  graphs <- names(seurat_obj@graphs)
  if (length(graphs) > 0) {
    for (graph in graphs) {
      components$obsp[[graph]] <- seurat_obj@graphs[[graph]]
    }
  }

  # Get additional assays as layers
  components$layers <- list()
  assays <- Seurat::Assays(seurat_obj)
  other_assays <- setdiff(assays, default_assay)
  if (length(other_assays) > 0) {
    for (assay in other_assays) {
      components$layers[[assay]] <- Matrix::t(Seurat::GetAssayData(seurat_obj, slot = "counts", assay = assay))
    }
  }

  # Get raw data from misc if present
  if (!is.null(seurat_obj@misc$raw)) {
    components$raw <- seurat_obj@misc$raw
  }

  # Get uns from misc
  if (!is.null(seurat_obj@misc$uns)) {
    components$uns <- seurat_obj@misc$uns
  } else {
    components$uns <- list()
  }

  # No varm or varp in standard Seurat
  components$varm <- list()
  components$varp <- list()

  return(components)
}
