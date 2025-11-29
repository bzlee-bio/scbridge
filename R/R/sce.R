# SingleCellExperiment object converters
# Part of the scio R package


#' Internal: Optimized transpose for sparse matrices
#'
#' Uses MatrixExtra::t_shallow() for O(1) transpose when available,
#' falls back to Matrix::t() otherwise.
#'
#' @param mat Sparse matrix (dgCMatrix)
#' @return Transposed matrix (dgRMatrix if MatrixExtra available, dgCMatrix otherwise)
.fast_transpose <- function(mat) {
  if (requireNamespace("MatrixExtra", quietly = TRUE)) {
    # O(1) transpose - returns dgRMatrix (CSR) which is semantically transposed
    # without copying data
    return(MatrixExtra::t_shallow(mat))
  } else {
    # Fallback to regular transpose
    return(Matrix::t(mat))
  }
}


#' Convert components list to SingleCellExperiment object
#'
#' @param components Named list with X, obs, var, obsm, obsp, etc.
#' @return SingleCellExperiment object
#' @export
create_sce_from_components <- function(components) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment package required. Install with: BiocManager::install('SingleCellExperiment')")
  }

  # Create assays list - transpose cells×genes to genes×cells
  # Use fast transpose (O(1) with MatrixExtra, regular otherwise)
  assays_list <- list(counts = .fast_transpose(components$X))

  # Add layers as additional assays
  if (length(components$layers) > 0) {
    for (layer_name in names(components$layers)) {
      assays_list[[layer_name]] <- .fast_transpose(components$layers[[layer_name]])
    }
  }

  # Add raw counts if present
  if (!is.null(components$raw)) {
    assays_list$raw_counts <- .fast_transpose(components$raw$X)
  }

  # Create SCE object
  sce_obj <- SingleCellExperiment::SingleCellExperiment(
    assays = assays_list,
    colData = components$obs,
    rowData = components$var
  )

  # Add embeddings (reducedDims)
  if (length(components$obsm) > 0) {
    for (emb_name in names(components$obsm)) {
      SingleCellExperiment::reducedDim(sce_obj, emb_name) <- components$obsm[[emb_name]]
    }
  }

  # Add cell-cell graphs to metadata
  if (length(components$obsp) > 0) {
    S4Vectors::metadata(sce_obj)$obsp <- components$obsp
  }

  # Add gene embeddings to metadata
  if (length(components$varm) > 0) {
    S4Vectors::metadata(sce_obj)$varm <- components$varm
  }

  # Add gene-gene graphs to metadata
  if (length(components$varp) > 0) {
    S4Vectors::metadata(sce_obj)$varp <- components$varp
  }

  # Add raw var to metadata if present
  if (!is.null(components$raw)) {
    S4Vectors::metadata(sce_obj)$raw_var <- components$raw$var
  }

  # Add uns metadata
  if (length(components$uns) > 0) {
    S4Vectors::metadata(sce_obj)$uns <- components$uns
  }

  return(sce_obj)
}


#' Extract components from SingleCellExperiment object
#'
#' @param sce_obj SingleCellExperiment object
#' @return Named list with all components
#' @export
extract_components_from_sce <- function(sce_obj) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment package required.")
  }

  components <- list()

  # Get count matrix (transpose from genes × cells to cells × genes)
  # Handle DelayedArray and other S4 matrix types by ensuring it's a proper matrix
  counts_matrix <- SummarizedExperiment::assay(sce_obj, "counts")
  if (inherits(counts_matrix, "DelayedArray")) {
    # For DelayedArray, realize to sparse or dense matrix first
    counts_matrix <- as(counts_matrix, "sparseMatrix")
  }
  # Use fast transpose (O(1) with MatrixExtra, or regular Matrix::t)
  components$X <- .fast_transpose(counts_matrix)

  # Get cell metadata
  components$obs <- as.data.frame(SummarizedExperiment::colData(sce_obj))

  # Get gene metadata
  components$var <- as.data.frame(SummarizedExperiment::rowData(sce_obj))

  # Get embeddings (reducedDims)
  components$obsm <- list()
  red_dim_names <- SingleCellExperiment::reducedDimNames(sce_obj)
  if (length(red_dim_names) > 0) {
    for (red_dim in red_dim_names) {
      components$obsm[[red_dim]] <- SingleCellExperiment::reducedDim(sce_obj, red_dim)
    }
  }

  # Get additional assays as layers
  components$layers <- list()
  assay_names <- SummarizedExperiment::assayNames(sce_obj)
  other_assays <- setdiff(assay_names, c("counts", "raw_counts"))
  if (length(other_assays) > 0) {
    for (assay in other_assays) {
      layer_matrix <- SummarizedExperiment::assay(sce_obj, assay)
      if (inherits(layer_matrix, "DelayedArray")) {
        layer_matrix <- as(layer_matrix, "sparseMatrix")
      }
      components$layers[[assay]] <- .fast_transpose(layer_matrix)
    }
  }

  # Get raw counts if present
  if ("raw_counts" %in% assay_names) {
    components$raw <- list()
    raw_matrix <- SummarizedExperiment::assay(sce_obj, "raw_counts")
    if (inherits(raw_matrix, "DelayedArray")) {
      raw_matrix <- as(raw_matrix, "sparseMatrix")
    }
    components$raw$X <- .fast_transpose(raw_matrix)

    # Get raw var from metadata if present
    if (!is.null(S4Vectors::metadata(sce_obj)$raw_var)) {
      components$raw$var <- S4Vectors::metadata(sce_obj)$raw_var
    } else {
      components$raw$var <- data.frame(row.names = colnames(components$raw$X))
    }
  }

  # Get obsp from metadata
  if (!is.null(S4Vectors::metadata(sce_obj)$obsp)) {
    components$obsp <- S4Vectors::metadata(sce_obj)$obsp
  } else {
    components$obsp <- list()
  }

  # Get varm from metadata
  if (!is.null(S4Vectors::metadata(sce_obj)$varm)) {
    components$varm <- S4Vectors::metadata(sce_obj)$varm
  } else {
    components$varm <- list()
  }

  # Get varp from metadata
  if (!is.null(S4Vectors::metadata(sce_obj)$varp)) {
    components$varp <- S4Vectors::metadata(sce_obj)$varp
  } else {
    components$varp <- list()
  }

  # Get uns from metadata
  if (!is.null(S4Vectors::metadata(sce_obj)$uns)) {
    components$uns <- S4Vectors::metadata(sce_obj)$uns
  } else {
    components$uns <- list()
  }

  return(components)
}
