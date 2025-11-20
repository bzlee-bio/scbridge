"""
Writers module - Save all 10 AnnData components to folder structure
Adapted from save_r_compatible_complete.py
"""
import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse as sp
from pathlib import Path
from typing import Optional, Dict
import warnings

from .formats import (
    save_mtx, save_parquet, save_json, save_tsv_gz,
    numpy_to_json_serializable
)


def save_to_folder(adata: ad.AnnData,
                   folder_path: Path,
                   compress: bool = True) -> Dict[str, Path]:
    """
    Save all AnnData components to a folder structure

    This saves ALL 10 components of AnnData:
    1. X - expression matrix (MTX)
    2. obs - cell metadata (Parquet)
    3. var - gene metadata (Parquet)
    4. obsm - cell embeddings (Parquet)
    5. varm - gene embeddings (Parquet)
    6. obsp - cell-cell graphs (MTX)
    7. varp - gene-gene graphs (MTX)
    8. layers - additional matrices (MTX)
    9. uns - unstructured metadata (JSON)
    10. raw - raw data (MTX + Parquet)

    Parameters:
    -----------
    adata : AnnData
        AnnData object to save
    folder_path : Path
        Output folder path
    compress : bool
        Whether to compress MTX files (default: True)

    Returns:
    --------
    dict : Dictionary mapping component names to saved file paths
    """
    folder_path = Path(folder_path)
    folder_path.mkdir(parents=True, exist_ok=True)

    saved_files = {}

    # =========================================================================
    # 1. X - Expression matrix (cells × genes, transposed to genes × cells for R)
    # =========================================================================
    # CRITICAL FIX: Don't use gzip on MTX - it's extremely slow (10-100x slower)
    # The tar archive itself provides sufficient compression
    matrix_file = "matrix.mtx"

    # Transpose for R compatibility (genes × cells)
    if sp.issparse(adata.X):
        X_transposed = adata.X.T
    else:
        X_transposed = sp.csr_matrix(adata.X.T)

    save_mtx(X_transposed, folder_path / matrix_file, compress=False)  # NO gzip - too slow!
    saved_files['X'] = matrix_file

    # =========================================================================
    # 2. obs - Cell metadata (barcodes + metadata)
    # =========================================================================
    # Save cell IDs as barcodes
    barcodes_file = "barcodes.tsv.gz" if compress else "barcodes.tsv"
    barcodes_df = pd.DataFrame(index=adata.obs_names)
    save_tsv_gz(barcodes_df, folder_path / barcodes_file, header=False)
    saved_files['barcodes'] = barcodes_file

    # Save cell metadata as Parquet
    save_parquet(adata.obs, folder_path / "obs.parquet")
    saved_files['obs'] = "obs.parquet"

    # =========================================================================
    # 3. var - Gene metadata (features + metadata)
    # =========================================================================
    # Save gene IDs as features
    features_file = "features.tsv.gz" if compress else "features.tsv"
    features_df = pd.DataFrame({
        'gene_id': adata.var_names,
        'gene_name': adata.var_names,
        'feature_type': 'Gene Expression'
    })
    features_df.to_csv(
        folder_path / features_file,
        sep='\t',
        header=False,
        index=False,
        compression='gzip' if compress else None
    )
    saved_files['features'] = features_file

    # Save gene metadata as Parquet
    if adata.var.shape[1] > 0:
        save_parquet(adata.var, folder_path / "var.parquet")
        saved_files['var'] = "var.parquet"

    # =========================================================================
    # 4. obsm - Cell embeddings (PCA, UMAP, etc.)
    # =========================================================================
    if len(adata.obsm.keys()) > 0:
        obsm_dir = folder_path / "obsm"
        obsm_dir.mkdir(exist_ok=True)

        for key in adata.obsm.keys():
            emb_array = adata.obsm[key]

            # Create DataFrame with proper column names
            emb_df = pd.DataFrame(
                emb_array,
                index=adata.obs_names,
                columns=[f"{key}_{i+1}" for i in range(emb_array.shape[1])]
            )

            save_parquet(emb_df, obsm_dir / f"{key}.parquet")
            saved_files[f'obsm_{key}'] = f"obsm/{key}.parquet"

    # =========================================================================
    # 5. varm - Gene embeddings
    # =========================================================================
    if len(adata.varm.keys()) > 0:
        varm_dir = folder_path / "varm"
        varm_dir.mkdir(exist_ok=True)

        for key in adata.varm.keys():
            varm_array = adata.varm[key]

            varm_df = pd.DataFrame(
                varm_array,
                index=adata.var_names,
                columns=[f"{key}_{i+1}" for i in range(varm_array.shape[1])]
            )

            save_parquet(varm_df, varm_dir / f"{key}.parquet")
            saved_files[f'varm_{key}'] = f"varm/{key}.parquet"

    # =========================================================================
    # 6. obsp - Cell-cell graphs (neighbors, distances)
    # =========================================================================
    if len(adata.obsp.keys()) > 0:
        obsp_dir = folder_path / "obsp"
        obsp_dir.mkdir(exist_ok=True)

        for key in adata.obsp.keys():
            graph_file = f"{key}.mtx"  # No gzip - too slow
            save_mtx(adata.obsp[key], obsp_dir / graph_file, compress=False)
            saved_files[f'obsp_{key}'] = f"obsp/{graph_file}"

    # =========================================================================
    # 7. varp - Gene-gene graphs
    # =========================================================================
    if len(adata.varp.keys()) > 0:
        varp_dir = folder_path / "varp"
        varp_dir.mkdir(exist_ok=True)

        for key in adata.varp.keys():
            graph_file = f"{key}.mtx"  # No gzip - too slow
            save_mtx(adata.varp[key], varp_dir / graph_file, compress=False)
            saved_files[f'varp_{key}'] = f"varp/{graph_file}"

    # =========================================================================
    # 8. layers - Additional matrices (raw counts, normalized, etc.)
    # =========================================================================
    if len(adata.layers.keys()) > 0:
        layers_dir = folder_path / "layers"
        layers_dir.mkdir(exist_ok=True)

        for key in adata.layers.keys():
            layer_file = f"{key}.mtx"  # No gzip - too slow

            # Transpose for R (genes × cells)
            if sp.issparse(adata.layers[key]):
                layer_transposed = adata.layers[key].T
            else:
                layer_transposed = sp.csr_matrix(adata.layers[key].T)

            save_mtx(layer_transposed, layers_dir / layer_file, compress=False)
            saved_files[f'layer_{key}'] = f"layers/{layer_file}"

    # =========================================================================
    # 9. raw - Raw data (if present)
    # =========================================================================
    if adata.raw is not None:
        raw_dir = folder_path / "raw"
        raw_dir.mkdir(exist_ok=True)

        # Save raw X matrix
        raw_matrix_file = "matrix.mtx"  # No gzip - too slow!

        if sp.issparse(adata.raw.X):
            raw_X_transposed = adata.raw.X.T
        else:
            raw_X_transposed = sp.csr_matrix(adata.raw.X.T)

        save_mtx(raw_X_transposed, raw_dir / raw_matrix_file, compress=compress)
        saved_files['raw_X'] = f"raw/{raw_matrix_file}"

        # Save raw var (gene metadata for raw counts)
        if adata.raw.var.shape[1] > 0:
            save_parquet(adata.raw.var, raw_dir / "var.parquet")
            saved_files['raw_var'] = "raw/var.parquet"

        # Save raw gene IDs
        raw_features_file = "features.tsv.gz" if compress else "features.tsv"
        raw_features_df = pd.DataFrame({
            'gene_id': adata.raw.var_names,
            'gene_name': adata.raw.var_names,
            'feature_type': 'Gene Expression'
        })
        raw_features_df.to_csv(
            raw_dir / raw_features_file,
            sep='\t',
            header=False,
            index=False,
            compression='gzip' if compress else None
        )
        saved_files['raw_features'] = f"raw/{raw_features_file}"

    # =========================================================================
    # 10. uns - Unstructured metadata
    # =========================================================================
    if len(adata.uns.keys()) > 0:
        uns_json_serializable = {}

        for key, value in adata.uns.items():
            try:
                # Convert numpy types to JSON-serializable
                serializable_value = numpy_to_json_serializable(value)
                # Test if it's JSON serializable
                import json
                json.dumps(serializable_value)
                uns_json_serializable[key] = serializable_value
            except (TypeError, ValueError):
                # Skip non-serializable objects
                warnings.warn(f"Skipping non-serializable uns entry: {key}")
                continue

        if uns_json_serializable:
            save_json(uns_json_serializable, folder_path / "uns.json")
            saved_files['uns'] = "uns.json"

    # =========================================================================
    # Create manifest file
    # =========================================================================
    manifest = {
        'format': 'scBridge v1.0',
        'created_by': 'scio.write',
        'dimensions': {
            'n_obs': int(adata.n_obs),
            'n_vars': int(adata.n_vars),
        },
        'components': {
            'X': True,
            'obs': True,
            'var': adata.var.shape[1] > 0,
            'obsm': list(adata.obsm.keys()),
            'varm': list(adata.varm.keys()),
            'obsp': list(adata.obsp.keys()),
            'varp': list(adata.varp.keys()),
            'layers': list(adata.layers.keys()),
            'raw': adata.raw is not None,
            'uns': list(adata.uns.keys()) if len(adata.uns.keys()) > 0 else [],
        },
        'files': saved_files
    }

    save_json(manifest, folder_path / "manifest.json")

    return saved_files
