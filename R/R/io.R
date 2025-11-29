# Main I/O module - Simple API for writing and reading data in .scio format
# Part of the scio R package

#' Write single-cell data to .scio folder
#'
#' Saves Seurat or SingleCellExperiment objects to .scio folder format
#' using binary sparse format for sparse matrices.
#' Works seamlessly with Python (loads as AnnData).
#'
#' v0.1.3: Supports both CSR and CSC sparse matrix formats (.npy files).
#' Data is stored in cells x genes orientation.
#'
#' @param object Seurat or SingleCellExperiment object
#' @param path Output .scio folder path (e.g., "data.scio")
#' @param sparse_format Format for sparse matrices: "csr" (default, fastest write)
#'   or "csc" (faster R read). Default is "csr" since scanpy outputs CSR.
#' @param overwrite Whether to overwrite existing folder (default: FALSE)
#' @param update Whether to perform incremental update using hash-based
#'   change detection (default: FALSE). Only changed components will be
#'   rewritten. If the existing file has no hash info, a full save is
#'   performed with hashes computed for future updates.
#'
#' @examples
#' \dontrun{
#' library(scio)
#'
#' # Save Seurat object (full write, default CSR format)
#' scio_write(seurat_obj, "data.scio")
#'
#' # Use CSC format for faster R reading
#' scio_write(seurat_obj, "data.scio", sparse_format = "csc")
#'
#' # Overwrite existing file
#' scio_write(seurat_obj, "data.scio", overwrite = TRUE)
#'
#' # Incremental update (only changed components)
#' scio_write(seurat_obj, "data.scio", update = TRUE)
#' }
#'
#' @export
scio_write <- function(object, path, sparse_format = "csr", overwrite = FALSE, update = FALSE) {
  # Validate sparse_format
  sparse_format <- match.arg(sparse_format, c("csr", "csc"))

  # Validate conflicting options
  if (overwrite && update) {
    stop("Cannot use both overwrite = TRUE and update = TRUE. Choose one.")
  }

  # Detect object type and extract components
  if (inherits(object, "Seurat")) {
    components <- extract_components_from_seurat(object)
  } else if (inherits(object, "SingleCellExperiment")) {
    components <- extract_components_from_sce(object)
  } else {
    stop("Unsupported object type. Must be Seurat or SingleCellExperiment.")
  }

  # Handle update mode
  if (update && dir.exists(path)) {
    # Load existing manifest
    manifest_path <- file.path(path, "manifest.json")
    if (!file.exists(manifest_path)) {
      stop(paste0("Invalid scio folder: manifest.json not found in ", path))
    }

    manifest <- load_json(manifest_path)

    if (is.null(manifest$hashes)) {
      # No hash info - do full save with hashes
      message("  No hash info found in existing file. Doing full save with hash computation...")
      save_to_folder(components, path, sparse_format = sparse_format, compute_hashes = TRUE)
      message("  \u2713 Data saved with hashes to ", basename(path))
    } else {
      # Incremental update - only write changed components
      message("  Performing incremental update...")
      update_folder(components, path, manifest, sparse_format = sparse_format)
      message("  \u2713 Incremental update complete for ", basename(path))
    }

    return(invisible(path))
  }

  # Check if folder exists for non-update mode
  if (dir.exists(path) && !overwrite && !update) {
    stop(paste0("Folder already exists: ", path, "\nUse overwrite = TRUE to replace it, or update = TRUE for incremental update."))
  }

  # Remove existing folder if overwrite is TRUE
  if (dir.exists(path) && overwrite) {
    unlink(path, recursive = TRUE)
  }

  # Save directly to .scio folder (skip hash computation for faster write)
  message("  Saving data to ", basename(path), "...")
  save_to_folder(components, path, sparse_format = sparse_format, compute_hashes = FALSE)
  message("  \u2713 Data saved to ", basename(path))

  invisible(path)
}


#' Read single-cell data from .scio folder
#'
#' Loads data from .scio folder and converts to Seurat or SingleCellExperiment.
#'
#' @param path Path to .scio folder (also supports legacy tar archives)
#' @param output Output format: "Seurat" or "SCE" (default: "Seurat")
#'
#' @return Seurat or SingleCellExperiment object
#'
#' @examples
#' \dontrun{
#' library(scio)
#'
#' # Load as Seurat
#' seurat_obj <- scio_read("data.scio", output = "Seurat")
#'
#' # Load as SingleCellExperiment
#' sce_obj <- scio_read("data.scio", output = "SCE")
#' }
#'
#' @export
scio_read <- function(path, output = "Seurat") {
  if (!file.exists(path)) {
    stop(paste0("File not found: ", path))
  }

  # Validate output format
  output <- match.arg(output, c("Seurat", "SCE"))

  # Check if it's a folder or a tar archive
  if (dir.exists(path)) {
    # It's a folder - load directly
    message("  Loading data from ", basename(path), "...")
    components <- load_from_folder(path)
    message("  \u2713 Data loaded")
  } else {
    # It's a file - try tar archive (for backward compatibility)
    message("  Extracting legacy tar archive...")
    tmpdir <- tempdir()

    # Extract tar archive
    untar(path, exdir = tmpdir)

    # Find the extracted folder
    folder_name <- tools::file_path_sans_ext(basename(path))
    data_folder <- file.path(tmpdir, folder_name)

    if (!dir.exists(data_folder)) {
      stop(paste0("Extracted folder not found: ", data_folder))
    }

    # Load from folder structure
    message("  Loading data...")
    components <- load_from_folder(data_folder)
    message("  \u2713 Data loaded")

    # Clean up temporary folder
    unlink(data_folder, recursive = TRUE)
  }

  # Convert to requested format
  if (output == "Seurat") {
    result <- create_seurat_from_components(components)
  } else if (output == "SCE") {
    result <- create_sce_from_components(components)
  }

  return(result)
}


#' Write components list to .scio folder (advanced usage)
#'
#' @param components Named list with all components
#' @param path Output .scio folder path
#' @param overwrite Whether to overwrite existing folder
#' @export
scio_write_components <- function(components, path, overwrite = FALSE) {
  # Check if folder exists
  if (dir.exists(path) && !overwrite) {
    stop(paste0("Folder already exists: ", path, "\nUse overwrite = TRUE to replace it."))
  }

  # Remove if exists
  if (dir.exists(path) && overwrite) {
    unlink(path, recursive = TRUE)
  }

  # Save to folder structure (skip hash computation for faster write)
  save_to_folder(components, path, compute_hashes = FALSE)

  invisible(path)
}


#' Read components list from .scio folder (advanced usage)
#'
#' @param path Path to .scio folder
#' @return Named list with all components
#' @export
scio_read_components <- function(path) {
  if (!file.exists(path)) {
    stop(paste0("Path not found: ", path))
  }

  if (dir.exists(path)) {
    # It's a folder - load directly
    components <- load_from_folder(path)
  } else {
    # Legacy tar archive support
    tmpdir <- tempdir()
    untar(path, exdir = tmpdir)
    folder_name <- tools::file_path_sans_ext(basename(path))
    data_folder <- file.path(tmpdir, folder_name)

    if (!dir.exists(data_folder)) {
      stop(paste0("Extracted folder not found: ", data_folder))
    }

    components <- load_from_folder(data_folder)
    unlink(data_folder, recursive = TRUE)
  }

  return(components)
}
