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
    Write AnnData object to tar file

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

    Format inside tar:
    - MTX for sparse matrices (universal R/Python compatibility)
    - Parquet for DataFrames (preserves dtypes, efficient for large data)
    - JSON for metadata

    Parameters:
    -----------
    adata : AnnData
        AnnData object to save
    path : str or Path
        Output tar file path (e.g., "data.tar")
    overwrite : bool
        Whether to overwrite existing file (default: False)
    compress : bool
        Whether to compress MTX files inside tar (default: True)

    Example:
    --------
    >>> import anndata as ad
    >>> import scbridge as sb
    >>> adata = ad.read_h5ad("data.h5ad")
    >>> sb.write(adata, "data.tar")
    """
    path = Path(path)

    # Check if file exists
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"File already exists: {path}\n"
            f"Use overwrite=True to replace it."
        )

    # Create temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        folder_name = path.stem  # Use tar filename as folder name
        data_folder = tmpdir / folder_name

        # Save to folder structure
        save_to_folder(adata, data_folder, compress=compress)

        # Create tar archive (uncompressed for speed)
        with tarfile.open(path, 'w') as tar:
            # Add the entire folder with arcname to preserve structure
            tar.add(data_folder, arcname=folder_name)


def read(path: Union[str, Path]) -> ad.AnnData:
    """
    Read AnnData object from tar file

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
        Path to tar file

    Returns:
    --------
    AnnData : Reconstructed AnnData object

    Example:
    --------
    >>> import scbridge as sb
    >>> adata = sb.read("data.tar")
    >>> print(adata)
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    # Create temporary directory for extraction
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Extract tar archive
        with tarfile.open(path, 'r') as tar:
            tar.extractall(tmpdir)

        # Find the extracted folder (should be only one)
        extracted_folders = [d for d in tmpdir.iterdir() if d.is_dir()]

        if len(extracted_folders) == 0:
            raise ValueError(f"No folder found in tar file: {path}")
        elif len(extracted_folders) > 1:
            raise ValueError(f"Multiple folders found in tar file: {path}")

        data_folder = extracted_folders[0]

        # Load from folder structure
        adata = load_from_folder(data_folder)

    return adata
