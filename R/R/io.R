# Main I/O module - Simple API for writing and reading data in .scio format
# Part of the scio R package

#' Write single-cell data to .scio folder
#'
#' Saves Seurat or SingleCellExperiment objects to .scio folder format.
#' Works seamlessly with Python (loads as AnnData).
#'
#' @param object Seurat or SingleCellExperiment object
#' @param path Output .scio folder path (e.g., "data.scio")
#' @param overwrite Whether to overwrite existing folder (default: FALSE)
#' @param compress Whether to compress MTX files (default: TRUE)
#'
#' @examples
#' \dontrun{
#' library(scio)
#'
#' # Save Seurat object
#' scio_write(seurat_obj, "data.scio")
#'
#' # Save SingleCellExperiment object
#' scio_write(sce_obj, "data.scio")
#' }
#'
#' @export
scio_write <- function(object, path, overwrite = FALSE, compress = TRUE) {
  # Check if folder exists
  if (dir.exists(path) && !overwrite) {
    stop(paste0("Folder already exists: ", path, "\nUse overwrite = TRUE to replace it."))
  }

  # Detect object type and extract components
  if (inherits(object, "Seurat")) {
    components <- extract_components_from_seurat(object)
  } else if (inherits(object, "SingleCellExperiment")) {
    components <- extract_components_from_sce(object)
  } else {
    stop("Unsupported object type. Must be Seurat or SingleCellExperiment.")
  }

  # Remove existing folder if overwrite is TRUE
  if (dir.exists(path) && overwrite) {
    unlink(path, recursive = TRUE)
  }

  # Save directly to .scio folder
  message("  Saving data to ", basename(path), "...")
  save_to_folder(components, path, compress = compress)
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


#' Write components list to .scio file (advanced usage)
#'
#' @param components Named list with all components
#' @param path Output .scio file path
#' @param overwrite Whether to overwrite existing file
#' @param compress Whether to compress MTX files
#' @export
scio_write_components <- function(components, path, overwrite = FALSE, compress = TRUE) {
  # Check if file exists
  if (file.exists(path) && !overwrite) {
    stop(paste0("File already exists: ", path, "\nUse overwrite = TRUE to replace it."))
  }

  # Create temporary directory
  tmpdir <- tempdir()
  folder_name <- tools::file_path_sans_ext(basename(path))
  data_folder <- file.path(tmpdir, folder_name)

  # Remove if exists
  if (dir.exists(data_folder)) {
    unlink(data_folder, recursive = TRUE)
  }

  # Save to folder structure
  save_to_folder(components, data_folder, compress = compress)

  # Create tar archive
  tar(path, files = folder_name, compression = "none", tar = "internal")

  # Clean up
  unlink(data_folder, recursive = TRUE)

  invisible(path)
}


#' Read components list from .scio file (advanced usage)
#'
#' @param path Path to .scio file
#' @return Named list with all components
#' @export
scio_read_components <- function(path) {
  if (!file.exists(path)) {
    stop(paste0("File not found: ", path))
  }

  # Create temporary directory for extraction
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
  components <- load_from_folder(data_folder)

  # Clean up
  unlink(data_folder, recursive = TRUE)

  return(components)
}
