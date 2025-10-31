"""
scbridge - Cross-platform single-cell RNA-seq data storage

A simple, efficient package for saving and loading AnnData objects in a format
that works seamlessly with both Python (AnnData) and R (Seurat/SingleCellExperiment).

Main functions:
- write(): Save AnnData to .scb file
- read(): Load AnnData from .scb file

Example:
    >>> import scbridge as sb
    >>> import anndata as ad
    >>>
    >>> # Save
    >>> adata = ad.read_h5ad("data.h5ad")
    >>> sb.write(adata, "data.scb")
    >>>
    >>> # Load
    >>> adata = sb.read("data.scb")
"""

from ._version import __version__
from .io import write, read

__all__ = ['write', 'read', '__version__']
