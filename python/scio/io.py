"""
Main I/O module - Simple API for writing and reading AnnData in tar format
"""
import anndata as ad
import tarfile
import tempfile
import shutil
from pathlib import Path
from typing import Union

from .writers import save_to_folder, update_folder
from .readers import load_from_folder
from .formats import load_json
from .lazy import LazyAnnData


def write(adata: ad.AnnData,
          path: Union[str, Path],
          sparse_format: str = 'csr',
          overwrite: bool = False,
          update: bool = False) -> None:
    """
    Write AnnData object to .scio folder

    Saves ALL 10 components of AnnData:
    - X (expression matrix)
    - obs (cell metadata)
    - var (gene metadata)
    - obsm (cell embeddings: PCA, UMAP, etc.)
    - varm (gene embeddings)
    - obsp (cell-cell graphs)
    - varp (gene-gene graphs)
    - layers (additional matrices)
    - raw (raw counts)
    - uns (unstructured metadata)

    Format inside .scio (v0.1.3):
    - Binary CSC/CSR (.npy) for sparse matrices (fast, cross-platform)
    - Parquet for DataFrames (preserves dtypes, efficient for large data)
    - JSON for metadata

    Parameters:
    -----------
    adata : AnnData
        AnnData object to save
    path : str or Path
        Output .scio folder path (e.g., "data.scio")
    sparse_format : str
        Format for sparse matrices: 'csr' (default, fastest write from scanpy)
        or 'csc' (faster R read). Default is 'csr' since scanpy outputs CSR.
    overwrite : bool
        Whether to overwrite existing folder (default: False)
    update : bool
        Whether to perform incremental update using hash-based change detection
        (default: False). Only changed components will be rewritten.

    Example:
    --------
    >>> import anndata as ad
    >>> import scio
    >>> adata = ad.read_h5ad("data.h5ad")
    >>> scio.write(adata, "data.scio")
    >>>
    >>> # Use CSC format for faster R reading
    >>> scio.write(adata, "data.scio", sparse_format='csc')
    >>>
    >>> # Incremental update (only changed components)
    >>> scio.write(adata, "data.scio", update=True)
    """
    path = Path(path)

    # Validate sparse_format
    if sparse_format not in ('csr', 'csc'):
        raise ValueError(f"sparse_format must be 'csr' or 'csc', got '{sparse_format}'")

    # Validate extension
    if not str(path).endswith('.scio'):
        raise ValueError(
            f"Invalid file extension. Only '.scio' is supported.\n"
            f"Got: {path.suffix}\n"
            f"Example: scio.write(adata, 'data.scio')"
        )

    # Validate conflicting options
    if overwrite and update:
        raise ValueError("Cannot use both overwrite=True and update=True. Choose one.")

    # Handle update mode
    if update and path.exists():
        manifest_path = path / "manifest.json"
        if not manifest_path.exists():
            raise ValueError(f"Invalid scio folder: manifest.json not found in {path}")

        manifest = load_json(manifest_path)

        if not manifest.get('hashes'):
            # No hash info - do full save with hashes
            print("  No hash info found in existing file. Doing full save with hash computation...", flush=True)
            save_to_folder(adata, path, sparse_format=sparse_format, compute_hashes=True)
            print(f"  ✓ Data saved with hashes to {path.name}", flush=True)
        else:
            # Incremental update - only write changed components
            print("  Performing incremental update...", flush=True)
            update_folder(adata, path, manifest, sparse_format=sparse_format)
            print(f"  ✓ Incremental update complete for {path.name}", flush=True)

        return

    # Check if folder exists for non-update mode
    if path.exists() and not overwrite and not update:
        raise FileExistsError(
            f"Folder already exists: {path}\n"
            f"Use overwrite=True to replace it, or update=True for incremental update."
        )

    # Remove existing folder if overwrite is True
    if path.exists() and overwrite:
        shutil.rmtree(path)

    print(f"  Saving data to {path.name}...", flush=True)
    # Compute hashes if update=True (for future incremental updates)
    # Skip hash computation on normal writes for speed
    save_to_folder(adata, path, sparse_format=sparse_format, compute_hashes=update)
    print(f"  ✓ Data saved to {path.name}", flush=True)


def read(path: Union[str, Path], lazy: bool = False) -> Union[ad.AnnData, LazyAnnData]:
    """
    Read AnnData object from .scio folder (or legacy tar archive)

    Loads ALL 10 components of AnnData:
    - X (expression matrix)
    - obs (cell metadata)
    - var (gene metadata)
    - obsm (cell embeddings)
    - varm (gene embeddings)
    - obsp (cell-cell graphs)
    - varp (gene-gene graphs)
    - layers (additional matrices)
    - raw (raw counts)
    - uns (unstructured metadata)

    Parameters:
    -----------
    path : str or Path
        Path to .scio folder (also supports legacy tar archives)
    lazy : bool
        If True, use lazy loading (only load components when accessed).
        This significantly improves performance for large datasets.
        Default: False (load all data immediately)

    Returns:
    --------
    AnnData or LazyAnnData :
        Reconstructed AnnData object (lazy=False) or
        LazyAnnData object (lazy=True)

    Example:
    --------
    >>> import scio
    >>> # Standard loading - loads all data
    >>> adata = scio.read("data.scio")
    >>>
    >>> # Lazy loading - only loads structure initially
    >>> adata = scio.read("data.scio", lazy=True)
    >>> print(adata.shape)  # Fast - no data loaded
    >>> cell_types = adata.obs['cell_type']  # Loads obs on first access
    >>> expression = adata.X  # Loads X on first access
    >>>
    >>> # Convert lazy to full AnnData when needed
    >>> full_adata = adata.to_memory()
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    # Check if it's a folder or a tar archive
    if path.is_dir():
        # It's a folder - load directly (or lazy load)
        if lazy:
            print(f"  Initializing lazy loading from {path.name}...", flush=True)
            adata = LazyAnnData(path)
            print(f"  ✓ Structure loaded (data will load on access)", flush=True)
            return adata
        else:
            print(f"  Loading data from {path.name}...", flush=True)
            adata = load_from_folder(path)
            print(f"  ✓ Data loaded", flush=True)
            return adata

    # It's a file - try to open as tar archive (for backward compatibility)
    if lazy:
        raise ValueError(
            "Lazy loading is not supported for legacy tar archives.\n"
            "Please use the folder format (.scio) for lazy loading."
        )

    try:
        print(f"  Extracting legacy tar archive...", flush=True)
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Extract tar archive
            with tarfile.open(path, 'r') as tar:
                tar.extractall(tmpdir)

            # Find the extracted folder (should be only one)
            extracted_folders = [d for d in tmpdir.iterdir() if d.is_dir()]

            if len(extracted_folders) == 0:
                raise ValueError(f"No folder found in archive: {path}")
            elif len(extracted_folders) > 1:
                raise ValueError(f"Multiple folders found in archive: {path}")

            data_folder = extracted_folders[0]

            # Load from folder structure
            print(f"  Loading data...", flush=True)
            adata = load_from_folder(data_folder)
            print(f"  ✓ Data loaded", flush=True)

        return adata

    except tarfile.ReadError as e:
        raise ValueError(
            f"Path is not a valid .scio folder or tar archive: {path}\n"
            f"Expected either a .scio folder or a legacy tar archive.\n"
            f"Error: {str(e)}"
        )
