"""
scio - Cross-platform single-cell RNA-seq data storage

A simple, efficient package for saving and loading AnnData objects in a format
that works seamlessly with both Python (AnnData) and R (Seurat/SingleCellExperiment).

Main functions:
- write(): Save AnnData to .scio file
- read(): Load AnnData from .scio file

Example:
    >>> import scio
    >>> import anndata as ad
    >>>
    >>> # Save
    >>> adata = ad.read_h5ad("data.h5ad")
    >>> scio.write(adata, "data.scio")
    >>>
    >>> # Load
    >>> adata = scio.read("data.scio")
"""

from ._version import __version__
from .io import write, read

__all__ = ['write', 'read', '__version__']
