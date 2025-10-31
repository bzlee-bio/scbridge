"""
Readers module - Load all 10 AnnData components from folder structure
"""
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional
import warnings

from .formats import (
    load_mtx, load_parquet, load_json, load_tsv_gz
)


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
    # 1. Load X matrix (transposed from genes × cells to cells × genes)
    # =========================================================================
    matrix_file = folder_path / manifest['files']['X']
    X = load_mtx(matrix_file, transpose=True)  # Transpose back to cells × genes

    # =========================================================================
    # 2. Load cell IDs (obs_names)
    # =========================================================================
    barcodes_file = folder_path / manifest['files']['barcodes']
    barcodes_df = load_tsv_gz(barcodes_file, header=False)
    obs_names = barcodes_df.iloc[:, 0].values

    # =========================================================================
    # 3. Load gene IDs (var_names)
    # =========================================================================
    features_file = folder_path / manifest['files']['features']
    features_df = pd.read_csv(
        features_file,
        sep='\t',
        header=None,
        compression='gzip' if str(features_file).endswith('.gz') else None
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
    obsm_keys = manifest['components']['obsm']
    if obsm_keys:
        for key in obsm_keys:
            file_path = folder_path / manifest['files'][f'obsm_{key}']
            emb_df = load_parquet(file_path)
            adata.obsm[key] = emb_df.values

    # =========================================================================
    # 7. Load varm (gene embeddings)
    # =========================================================================
    varm_keys = manifest['components']['varm']
    if varm_keys:
        for key in varm_keys:
            file_path = folder_path / manifest['files'][f'varm_{key}']
            varm_df = load_parquet(file_path)
            adata.varm[key] = varm_df.values

    # =========================================================================
    # 8. Load obsp (cell-cell graphs)
    # =========================================================================
    obsp_keys = manifest['components']['obsp']
    if obsp_keys:
        for key in obsp_keys:
            file_path = folder_path / manifest['files'][f'obsp_{key}']
            adata.obsp[key] = load_mtx(file_path, transpose=False)

    # =========================================================================
    # 9. Load varp (gene-gene graphs)
    # =========================================================================
    varp_keys = manifest['components']['varp']
    if varp_keys:
        for key in varp_keys:
            file_path = folder_path / manifest['files'][f'varp_{key}']
            adata.varp[key] = load_mtx(file_path, transpose=False)

    # =========================================================================
    # 10. Load layers (additional matrices)
    # =========================================================================
    layer_keys = manifest['components']['layers']
    if layer_keys:
        for key in layer_keys:
            file_path = folder_path / manifest['files'][f'layer_{key}']
            # Transpose back to cells × genes
            adata.layers[key] = load_mtx(file_path, transpose=True)

    # =========================================================================
    # 11. Load raw data (if present)
    # =========================================================================
    if manifest['components']['raw']:
        raw_dir = folder_path / "raw"

        # Load raw X matrix
        raw_matrix_file = raw_dir / manifest['files']['raw_X'].split('/')[-1]
        raw_X = load_mtx(raw_matrix_file, transpose=True)  # Transpose to cells × genes

        # Load raw gene IDs
        raw_features_file = raw_dir / manifest['files']['raw_features'].split('/')[-1]
        raw_features_df = pd.read_csv(
            raw_features_file,
            sep='\t',
            header=None,
            compression='gzip' if str(raw_features_file).endswith('.gz') else None
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
