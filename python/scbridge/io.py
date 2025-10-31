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
    Write AnnData object to .scb file

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

    Format inside .scb:
    - MTX for sparse matrices (universal R/Python compatibility)
    - Parquet for DataFrames (preserves dtypes, efficient for large data)
    - JSON for metadata

    Parameters:
    -----------
    adata : AnnData
        AnnData object to save
    path : str or Path
        Output .scb file path (e.g., "data.scb")
    overwrite : bool
        Whether to overwrite existing file (default: False)
    compress : bool
        Whether to compress MTX files inside tar (default: True)

    Example:
    --------
    >>> import anndata as ad
    >>> import scbridge as sb
    >>> adata = ad.read_h5ad("data.h5ad")
    >>> sb.write(adata, "data.scb")
    """
    path = Path(path)

    # Validate extension
    if not str(path).endswith('.scb'):
        raise ValueError(
            f"Invalid file extension. Only '.scb' is supported.\n"
            f"Got: {path.suffix}\n"
            f"Example: sb.write(adata, 'data.scb')"
        )

    # Check if file exists
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"File already exists: {path}\n"
            f"Use overwrite=True to replace it."
        )

    # Create temporary folder for data
    folder_path = path.parent / path.stem

    print(f"  Saving data to folder...", flush=True)
    save_to_folder(adata, folder_path, compress=compress)
    print(f"  ✓ Data saved", flush=True)

    # Create uncompressed tar archive (fast - uses system tar command)
    print(f"  Creating .scb archive...", flush=True)
    import subprocess

    result = subprocess.run(
        ['tar', '-cf', str(path), '-C', str(folder_path.parent), folder_path.name],
        capture_output=True,
        text=True
    )

    if result.returncode == 0:
        # Remove folder after successful archive creation
        shutil.rmtree(folder_path)
        archive_size = path.stat().st_size / 1024 / 1024 / 1024  # GB
        print(f"  ✓ Archive created: {path.name} ({archive_size:.2f} GB)", flush=True)
    else:
        print(f"  ⚠ Archive creation failed: {result.stderr}", flush=True)
        print(f"  Keeping folder: {folder_path}", flush=True)


def read(path: Union[str, Path]) -> ad.AnnData:
    """
    Read AnnData object from .scb file (or any tar archive)

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
        Path to .scb file or any valid tar archive in scBridge format
        (supports .scb, .scbridge, .tar, or no extension)

    Returns:
    --------
    AnnData : Reconstructed AnnData object

    Example:
    --------
    >>> import scbridge as sb
    >>> adata = sb.read("data.scb")
    >>> # Also works with legacy formats:
    >>> adata = sb.read("data.scbridge")
    >>> adata = sb.read("data.tar")
    >>> print(adata)
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    # Try to open as tar archive - will raise an error if not a valid tar file
    try:
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
            adata = load_from_folder(data_folder)

        return adata

    except tarfile.ReadError as e:
        raise ValueError(
            f"File is not a valid tar archive: {path}\n"
            f"Expected a .scb file (tar archive format).\n"
            f"Error: {str(e)}"
        )
