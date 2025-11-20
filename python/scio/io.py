"""
Main I/O module - Simple API for writing and reading AnnData in tar format
"""
import anndata as ad
import tarfile
import tempfile
import shutil
from pathlib import Path
from typing import Union

from .writers import save_to_folder
from .readers import load_from_folder


def write(adata: ad.AnnData,
          path: Union[str, Path],
          overwrite: bool = False,
          compress: bool = True) -> None:
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

    Format inside .scio:
    - MTX for sparse matrices (universal R/Python compatibility)
    - Parquet for DataFrames (preserves dtypes, efficient for large data)
    - JSON for metadata

    Parameters:
    -----------
    adata : AnnData
        AnnData object to save
    path : str or Path
        Output .scio folder path (e.g., "data.scio")
    overwrite : bool
        Whether to overwrite existing folder (default: False)
    compress : bool
        Whether to compress MTX files (default: True)

    Example:
    --------
    >>> import anndata as ad
    >>> import scio as sb
    >>> adata = ad.read_h5ad("data.h5ad")
    >>> sb.write(adata, "data.scio")
    """
    path = Path(path)

    # Validate extension
    if not str(path).endswith('.scio'):
        raise ValueError(
            f"Invalid file extension. Only '.scio' is supported.\n"
            f"Got: {path.suffix}\n"
            f"Example: sb.write(adata, 'data.scio')"
        )

    # Check if folder exists
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Folder already exists: {path}\n"
            f"Use overwrite=True to replace it."
        )

    # Remove existing folder if overwrite is True
    if path.exists() and overwrite:
        shutil.rmtree(path)

    print(f"  Saving data to {path.name}...", flush=True)
    save_to_folder(adata, path, compress=compress)
    print(f"  ✓ Data saved to {path.name}", flush=True)


def read(path: Union[str, Path]) -> ad.AnnData:
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

    Returns:
    --------
    AnnData : Reconstructed AnnData object

    Example:
    --------
    >>> import scio as sb
    >>> adata = sb.read("data.scio")
    >>> # Also works with legacy tar archives:
    >>> adata = sb.read("data.scio")
    >>> print(adata)
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    # Check if it's a folder or a tar archive
    if path.is_dir():
        # It's a folder - load directly
        print(f"  Loading data from {path.name}...", flush=True)
        adata = load_from_folder(path)
        print(f"  ✓ Data loaded", flush=True)
        return adata

    # It's a file - try to open as tar archive (for backward compatibility)
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
