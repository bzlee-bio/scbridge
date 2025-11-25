# Format handlers for MTX, Parquet, and JSON files
# Part of the scio R package

#' Compute hash of an R object for change detection
#'
#' Uses MD5 hash for fast computation. Used for incremental updates.
#'
#' @param obj R object to hash
#' @return MD5 hash string
compute_hash <- function(obj) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("digest package required. Install with: install.packages('digest')")
  }
  digest::digest(obj, algo = "md5")
}

#' Load sparse matrix from MTX format
#'
#' @param mtx_file Path to MTX file (can be gzipped)
#' @return Sparse matrix
#' @importFrom Matrix readMM
load_mtx <- function(mtx_file) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package required. Install with: install.packages('Matrix')")
  }

  if (grepl("\\.gz$", mtx_file)) {
    con <- gzfile(mtx_file, "rb")
    mat <- Matrix::readMM(con)
    close(con)
  } else {
    mat <- Matrix::readMM(mtx_file)
  }

  return(mat)
}


#' Save sparse matrix to MTX format
#'
#' @param matrix Sparse matrix to save
#' @param file_path Output file path
#' @param compress Whether to gzip compress
#' @importFrom Matrix writeMM
save_mtx <- function(matrix, file_path, compress = TRUE) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package required. Install with: install.packages('Matrix')")
  }

  if (compress) {
    # Matrix::writeMM doesn't support connections, so write uncompressed first
    temp_file <- tempfile(fileext = ".mtx")
    Matrix::writeMM(matrix, temp_file)

    # Read and compress
    if (!grepl("\\.gz$", file_path)) {
      file_path <- paste0(file_path, ".gz")
    }

    con_in <- file(temp_file, "rb")
    con_out <- gzfile(file_path, "wb")
    writeBin(readBin(con_in, "raw", file.info(temp_file)$size), con_out)
    close(con_in)
    close(con_out)

    unlink(temp_file)
  } else {
    Matrix::writeMM(matrix, file_path)
  }

  return(file_path)
}


#' Load DataFrame from Parquet format
#'
#' @param parquet_file Path to Parquet file
#' @return Data frame
load_parquet <- function(parquet_file) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("arrow package required. Install with: install.packages('arrow')")
  }

  df <- arrow::read_parquet(parquet_file)
  return(as.data.frame(df))
}


#' Save DataFrame to Parquet format
#'
#' @param df Data frame to save
#' @param file_path Output file path
save_parquet <- function(df, file_path) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("arrow package required. Install with: install.packages('arrow')")
  }

  arrow::write_parquet(df, file_path)
  return(file_path)
}


#' Load list from JSON format
#'
#' @param json_file Path to JSON file
#' @return Named list
load_json <- function(json_file) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite package required. Install with: install.packages('jsonlite')")
  }

  data <- jsonlite::fromJSON(json_file)
  return(data)
}


#' Save list to JSON format
#'
#' @param data Named list to save
#' @param file_path Output file path
#' @param pretty Whether to pretty-print JSON
save_json <- function(data, file_path, pretty = TRUE) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite package required. Install with: install.packages('jsonlite')")
  }

  jsonlite::write_json(data, file_path, pretty = pretty, auto_unbox = TRUE)
  return(file_path)
}


#' Load sparse matrix from binary CSC format (v2.x)
#'
#' This loads the binary CSC format written by Python (v2.0+).
#' Files: {base_path}.data.npy, {base_path}.indices.npy, {base_path}.indptr.npy, {base_path}.shape.json
#'
#' Uses reticulate+numpy for reliable large array loading. Falls back to RcppCNPy if reticulate
#' is not available.
#'
#' @param base_path Base path to binary files (without extension)
#' @param transpose Whether to transpose after loading (default: FALSE)
#' @return dgCMatrix (or dgRMatrix if transpose=TRUE with MatrixExtra available)
#' @importFrom Matrix sparseMatrix
load_sparse_binary <- function(base_path, transpose = FALSE) {
  # Load shape metadata
  shape_file <- paste0(base_path, ".shape.json")
  if (!file.exists(shape_file)) {
    stop(paste("Shape file not found:", shape_file))
  }
  shape_meta <- load_json(shape_file)

  # Try reticulate+numpy first (most reliable for large arrays)
  if (requireNamespace("reticulate", quietly = TRUE)) {
    np <- tryCatch({
      reticulate::import("numpy", convert = FALSE)
    }, error = function(e) NULL)

    if (!is.null(np)) {
      # Load using numpy (handles int32 and large arrays reliably)
      data <- reticulate::py_to_r(np$load(paste0(base_path, ".data.npy")))
      indices <- reticulate::py_to_r(np$load(paste0(base_path, ".indices.npy")))
      indptr <- reticulate::py_to_r(np$load(paste0(base_path, ".indptr.npy")))
    } else {
      # Python/numpy not available, fall back to RcppCNPy
      data <- NULL
    }
  } else {
    data <- NULL
  }

  # Fallback to RcppCNPy if reticulate failed
  if (is.null(data)) {
    if (!requireNamespace("RcppCNPy", quietly = TRUE)) {
      stop("Either reticulate+numpy or RcppCNPy package required for binary format.\n",
           "Install with: install.packages('reticulate') or install.packages('RcppCNPy')")
    }
    data <- RcppCNPy::npyLoad(paste0(base_path, ".data.npy"))
    indices <- RcppCNPy::npyLoad(paste0(base_path, ".indices.npy"))
    indptr <- RcppCNPy::npyLoad(paste0(base_path, ".indptr.npy"))
  }

  # Reconstruct dgCMatrix (CSC format)
  # Note: R uses 0-based indices internally for dgCMatrix, same as Python
  mat <- new("dgCMatrix",
    x = as.numeric(data),
    i = as.integer(indices),
    p = as.integer(indptr),
    Dim = as.integer(shape_meta$shape)
  )

  if (transpose) {
    # Try to use MatrixExtra::t_shallow() for O(1) transpose if available
    # t_shallow() returns dgRMatrix (CSR) which is semantically transposed without data copy
    if (requireNamespace("MatrixExtra", quietly = TRUE)) {
      mat <- MatrixExtra::t_shallow(mat)
    } else {
      # Fallback to regular transpose (creates new matrix)
      mat <- Matrix::t(mat)
    }
  }

  return(mat)
}


#' Save sparse matrix to binary CSC format
#'
#' Saves a sparse matrix in binary CSC format for fast loading.
#' Creates: {base_path}.data.npy, {base_path}.indices.npy, {base_path}.indptr.npy, {base_path}.shape.json
#'
#' @param matrix Sparse matrix to save (will be converted to CSC if needed)
#' @param base_path Base output path (without extension)
#' @param orientation Data orientation metadata (e.g., "cells_x_genes" or "genes_x_cells")
#' @return Base path to saved files
save_sparse_binary <- function(matrix, base_path, orientation = NULL) {
  if (!requireNamespace("RcppCNPy", quietly = TRUE)) {
    stop("RcppCNPy package required for binary format. Install with: install.packages('RcppCNPy')")
  }

  # Convert to dgCMatrix (CSC format) if needed
  if (!inherits(matrix, "dgCMatrix")) {
    matrix <- as(matrix, "dgCMatrix")
  }

  # Save CSC components as binary numpy files
  RcppCNPy::npySave(paste0(base_path, ".data.npy"), matrix@x)
  RcppCNPy::npySave(paste0(base_path, ".indices.npy"), matrix@i)
  RcppCNPy::npySave(paste0(base_path, ".indptr.npy"), matrix@p)

  # Save shape metadata
  shape_meta <- list(
    shape = matrix@Dim,
    dtype = "float64",
    nnz = length(matrix@x),
    format = "csc_binary"
  )

  if (!is.null(orientation)) {
    shape_meta$orientation <- orientation
  }

  save_json(shape_meta, paste0(base_path, ".shape.json"))

  return(base_path)
}
