"""
Readers module - Load all 10 AnnData components from folder structure
"""
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, List, Union
import warnings

from .formats import (
    load_mtx, load_parquet, load_json, load_tsv_gz,
    load_sparse_parquet
)


def _ensure_list(value: Union[str, List, None]) -> List:
    """
    Ensure value is a list (handles R's auto_unbox converting single elements to strings)
    """
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, list):
        return value
    return list(value)


def _detect_gzip(file_path: Path) -> bool:
    """
    Detect if a file is gzipped by checking its magic bytes (0x1f 0x8b)

    Parameters:
    -----------
    file_path : Path
        Path to file to check

    Returns:
    --------
    bool : True if file is gzipped, False otherwise
    """
    try:
        with open(file_path, 'rb') as f:
            magic = f.read(2)
            return len(magic) == 2 and magic[0] == 0x1f and magic[1] == 0x8b
    except Exception:
        return False


def _load_sparse_matrix(folder_path: Path, file_key: str, manifest: dict,
                        transpose: bool = False, is_expression_matrix: bool = False):
    """
    Load a sparse matrix from either binary CSC (v2.x) or MTX (v1.x) format.

    Automatically detects the format based on file path:
    - v0.1.2+: base path (e.g., "matrix") -> loads matrix.data.npy, matrix.indices.npy, etc.
             with cells×genes orientation - no transpose needed for expression matrices
    - v2.0: same binary format but genes×cells orientation - need transpose
    - v1.x: mtx file path (e.g., "matrix.mtx.gz") -> loads via load_mtx

    Parameters:
    -----------
    is_expression_matrix : bool
        If True, this is X, layer, or raw.X matrix that should be cells×genes.
        For v0.1.2 cells_x_genes format, no transpose is needed.
        For v2.0 genes_x_cells format, transpose is needed.
    """
    file_base = manifest['files'][file_key]
    file_path = folder_path / file_base
    file_str = str(file_base)

    # Check orientation from manifest (v0.1.2+)
    orientation = manifest.get('orientation', 'genes_x_cells')  # Default for legacy

    # Determine if we need to transpose
    if is_expression_matrix:
        # For expression matrices (X, layers, raw.X)
        if orientation == 'cells_x_genes':
            # v0.1.2: already in cells×genes format - no transpose needed
            need_transpose = False
        else:
            # v2.0 or v1.x: genes×cells format - need transpose
            need_transpose = True
    else:
        # For non-expression matrices (obsp, varp), use the provided transpose flag
        need_transpose = transpose

    # v2.0+ binary CSC format: base path points to files like *.data.npy, *.shape.json
    shape_file = folder_path / (file_base + '.shape.json')
    if shape_file.exists():
        # v2.0+ binary CSC format
        return load_sparse_parquet(file_path, shape_file=shape_file, transpose=need_transpose)
    elif file_str.endswith('.mtx') or file_str.endswith('.mtx.gz'):
        # v1.x MTX format (always genes×cells, need transpose for expression matrices)
        return load_mtx(file_path, transpose=need_transpose)
    else:
        # Fallback: try v2.0+ binary format
        return load_sparse_parquet(file_path, transpose=need_transpose)


def load_from_folder(folder_path: Path) -> ad.AnnData:
    """
    Load all AnnData components from a folder structure

    This loads ALL 10 components of AnnData:
    1. X - expression matrix
    2. obs - cell metadata
    3. var - gene metadata
    4. obsm - cell embeddings
    5. varm - gene embeddings
    6. obsp - cell-cell graphs
    7. varp - gene-gene graphs
    8. layers - additional matrices
    9. uns - unstructured metadata
    10. raw - raw data

    Parameters:
    -----------
    folder_path : Path
        Folder containing saved data

    Returns:
    --------
    AnnData : Reconstructed AnnData object
    """
    folder_path = Path(folder_path)

    if not folder_path.exists():
        raise FileNotFoundError(f"Folder not found: {folder_path}")

    # Load manifest to understand structure
    manifest_path = folder_path / "manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Manifest file not found: {manifest_path}")

    manifest = load_json(manifest_path)

    # =========================================================================
    # 1. Load X matrix (cells × genes)
    # v0.1.2: already cells×genes, no transpose needed
    # v2.0/v1.x: genes×cells, transpose needed (handled by is_expression_matrix flag)
    # =========================================================================
    X = _load_sparse_matrix(folder_path, 'X', manifest, is_expression_matrix=True)

    # =========================================================================
    # 2. Load cell IDs (obs_names)
    # =========================================================================
    barcodes_file = folder_path / manifest['files']['barcodes']
    # Detect compression by checking file magic bytes (not just extension)
    compression = 'gzip' if _detect_gzip(barcodes_file) else None
    barcodes_df = pd.read_csv(barcodes_file, sep='\t', header=None, compression=compression)
    obs_names = barcodes_df.iloc[:, 0].values

    # =========================================================================
    # 3. Load gene IDs (var_names)
    # =========================================================================
    features_file = folder_path / manifest['files']['features']
    # Detect compression by checking file magic bytes (not just extension)
    compression = 'gzip' if _detect_gzip(features_file) else None
    features_df = pd.read_csv(
        features_file,
        sep='\t',
        header=None,
        compression=compression
    )
    var_names = features_df.iloc[:, 0].values

    # =========================================================================
    # 4. Load obs (cell metadata)
    # =========================================================================
    obs_file = folder_path / manifest['files']['obs']
    obs = load_parquet(obs_file)

    # =========================================================================
    # 5. Load var (gene metadata)
    # =========================================================================
    if 'var' in manifest['files']:
        var_file = folder_path / manifest['files']['var']
        var = load_parquet(var_file)
    else:
        # Create empty var with correct index
        var = pd.DataFrame(index=var_names)

    # =========================================================================
    # Create base AnnData object
    # =========================================================================
    adata = ad.AnnData(X=X, obs=obs, var=var)

    # Ensure correct index names
    adata.obs_names = obs_names
    adata.var_names = var_names

    # =========================================================================
    # 6. Load obsm (cell embeddings)
    # =========================================================================
    obsm_keys = _ensure_list(manifest['components']['obsm'])
    for key in obsm_keys:
        file_path = folder_path / manifest['files'][f'obsm_{key}']
        emb_df = load_parquet(file_path)
        adata.obsm[key] = emb_df.values

    # =========================================================================
    # 7. Load varm (gene embeddings)
    # =========================================================================
    varm_keys = _ensure_list(manifest['components']['varm'])
    for key in varm_keys:
        file_path = folder_path / manifest['files'][f'varm_{key}']
        varm_df = load_parquet(file_path)
        adata.varm[key] = varm_df.values

    # =========================================================================
    # 8. Load obsp (cell-cell graphs)
    # =========================================================================
    obsp_keys = _ensure_list(manifest['components']['obsp'])
    for key in obsp_keys:
        adata.obsp[key] = _load_sparse_matrix(folder_path, f'obsp_{key}', manifest, transpose=False)

    # =========================================================================
    # 9. Load varp (gene-gene graphs)
    # =========================================================================
    varp_keys = _ensure_list(manifest['components']['varp'])
    for key in varp_keys:
        adata.varp[key] = _load_sparse_matrix(folder_path, f'varp_{key}', manifest, transpose=False)

    # =========================================================================
    # 10. Load layers (additional matrices)
    # v0.1.2: already cells×genes, no transpose needed
    # v2.0/v1.x: genes×cells, transpose needed (handled by is_expression_matrix flag)
    # =========================================================================
    layer_keys = _ensure_list(manifest['components']['layers'])
    for key in layer_keys:
        adata.layers[key] = _load_sparse_matrix(folder_path, f'layer_{key}', manifest, is_expression_matrix=True)

    # =========================================================================
    # 11. Load raw data (if present)
    # v0.1.2: raw.X already cells×genes, no transpose needed
    # v2.0/v1.x: genes×cells, transpose needed (handled by is_expression_matrix flag)
    # =========================================================================
    if manifest['components']['raw']:
        raw_dir = folder_path / "raw"

        # Load raw X matrix
        raw_X = _load_sparse_matrix(folder_path, 'raw_X', manifest, is_expression_matrix=True)

        # Load raw gene IDs
        raw_features_file = raw_dir / manifest['files']['raw_features'].split('/')[-1]
        # Detect compression by checking file magic bytes (not just extension)
        compression = 'gzip' if _detect_gzip(raw_features_file) else None
        raw_features_df = pd.read_csv(
            raw_features_file,
            sep='\t',
            header=None,
            compression=compression
        )
        raw_var_names = raw_features_df.iloc[:, 0].values

        # Load raw var metadata if present
        if 'raw_var' in manifest['files']:
            raw_var_file = folder_path / manifest['files']['raw_var']
            raw_var = load_parquet(raw_var_file)
        else:
            raw_var = pd.DataFrame(index=raw_var_names)

        # Create raw AnnData
        raw_adata = ad.AnnData(X=raw_X, var=raw_var)
        raw_adata.var_names = raw_var_names
        raw_adata.obs_names = obs_names  # Same cells as main adata

        # Assign to main adata
        adata.raw = raw_adata

    # =========================================================================
    # 12. Load uns (unstructured metadata)
    # =========================================================================
    if 'uns' in manifest['files']:
        uns_file = folder_path / manifest['files']['uns']
        adata.uns = load_json(uns_file)

    return adata
