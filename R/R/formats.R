# Format handlers for MTX, Parquet, and JSON files
# Part of the scBridge R package

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
    mat <- Matrix::readMM(gzfile(mtx_file))
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

  if (compress && !grepl("\\.gz$", file_path)) {
    file_path <- paste0(file_path, ".gz")
  }

  if (compress) {
    con <- gzfile(file_path, "wb")
    Matrix::writeMM(matrix, con)
    close(con)
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
