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
    save_sparse_parquet,
    numpy_to_json_serializable, compute_hash, load_json
)


def save_to_folder(adata: ad.AnnData,
                   folder_path: Path,
                   compress: bool = True,
                   compute_hashes: bool = True) -> Dict[str, Path]:
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
    compute_hashes : bool
        Whether to compute and store MD5 hashes for incremental updates (default: True)

    Returns:
    --------
    dict : Dictionary mapping component names to saved file paths
    """
    folder_path = Path(folder_path)
    folder_path.mkdir(parents=True, exist_ok=True)

    saved_files = {}
    hashes = {}  # Store hashes for incremental updates

    # =========================================================================
    # 1. X - Expression matrix (cells × genes) - store as CSC directly
    # =========================================================================
    matrix_base = "matrix"

    # v0.1.2: Store in cells×genes orientation (Python's native format)
    # This eliminates expensive transpose operations on write
    # R will transpose on read using O(1) MatrixExtra::t_shallow()
    if sp.issparse(adata.X):
        if sp.isspmatrix_csc(adata.X):
            # Already CSC cells×genes - store directly
            X_to_store = adata.X
        elif sp.isspmatrix_csr(adata.X):
            # CSR.T gives CSC (cells×genes transposed = genes×cells, but shape stays same)
            # Actually we want cells×genes, so just convert CSR to CSC
            X_to_store = adata.X.tocsc()
        else:
            X_to_store = sp.csc_matrix(adata.X)
    else:
        X_to_store = sp.csc_matrix(adata.X)

    save_sparse_parquet(X_to_store, folder_path / matrix_base, orientation='cells_x_genes')
    saved_files['X'] = matrix_base
    if compute_hashes:
        hashes['X'] = compute_hash(adata.X)

    # =========================================================================
    # 2. obs - Cell metadata (barcodes + metadata)
    # =========================================================================
    # Save cell IDs as barcodes (always gzip - small file, fast compression)
    barcodes_file = "barcodes.tsv.gz"
    barcodes_df = pd.DataFrame(index=adata.obs_names)
    save_tsv_gz(barcodes_df, folder_path / barcodes_file, header=False)
    saved_files['barcodes'] = barcodes_file

    # Save cell metadata as Parquet
    save_parquet(adata.obs, folder_path / "obs.parquet")
    saved_files['obs'] = "obs.parquet"
    if compute_hashes:
        hashes['obs'] = compute_hash(adata.obs)

    # =========================================================================
    # 3. var - Gene metadata (features + metadata)
    # =========================================================================
    # Save gene IDs as features (always gzip - small file, fast compression)
    features_file = "features.tsv.gz"
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
        compression='gzip'
    )
    saved_files['features'] = features_file

    # Save gene metadata as Parquet
    if adata.var.shape[1] > 0:
        save_parquet(adata.var, folder_path / "var.parquet")
        saved_files['var'] = "var.parquet"
        if compute_hashes:
            hashes['var'] = compute_hash(adata.var)

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
            if compute_hashes:
                hashes[f'obsm_{key}'] = compute_hash(adata.obsm[key])

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
            if compute_hashes:
                hashes[f'varm_{key}'] = compute_hash(adata.varm[key])

    # =========================================================================
    # 6. obsp - Cell-cell graphs (neighbors, distances)
    # =========================================================================
    if len(adata.obsp.keys()) > 0:
        obsp_dir = folder_path / "obsp"
        obsp_dir.mkdir(exist_ok=True)

        for key in adata.obsp.keys():
            graph_base = key
            save_sparse_parquet(adata.obsp[key], obsp_dir / graph_base)
            saved_files[f'obsp_{key}'] = f"obsp/{graph_base}"
            if compute_hashes:
                hashes[f'obsp_{key}'] = compute_hash(adata.obsp[key])

    # =========================================================================
    # 7. varp - Gene-gene graphs
    # =========================================================================
    if len(adata.varp.keys()) > 0:
        varp_dir = folder_path / "varp"
        varp_dir.mkdir(exist_ok=True)

        for key in adata.varp.keys():
            graph_base = key
            save_sparse_parquet(adata.varp[key], varp_dir / graph_base)
            saved_files[f'varp_{key}'] = f"varp/{graph_base}"
            if compute_hashes:
                hashes[f'varp_{key}'] = compute_hash(adata.varp[key])

    # =========================================================================
    # 8. layers - Additional matrices (raw counts, normalized, etc.)
    # =========================================================================
    if len(adata.layers.keys()) > 0:
        layers_dir = folder_path / "layers"
        layers_dir.mkdir(exist_ok=True)

        for key in adata.layers.keys():
            layer_base = key

            # v0.1.2: Store in cells×genes orientation (no transpose)
            layer = adata.layers[key]
            if sp.issparse(layer):
                if sp.isspmatrix_csc(layer):
                    layer_to_store = layer
                elif sp.isspmatrix_csr(layer):
                    layer_to_store = layer.tocsc()
                else:
                    layer_to_store = sp.csc_matrix(layer)
            else:
                layer_to_store = sp.csc_matrix(layer)

            save_sparse_parquet(layer_to_store, layers_dir / layer_base, orientation='cells_x_genes')
            saved_files[f'layer_{key}'] = f"layers/{layer_base}"
            if compute_hashes:
                hashes[f'layer_{key}'] = compute_hash(adata.layers[key])

    # =========================================================================
    # 9. raw - Raw data (if present)
    # =========================================================================
    if adata.raw is not None:
        raw_dir = folder_path / "raw"
        raw_dir.mkdir(exist_ok=True)

        # Save raw X matrix - v0.1.2: cells×genes orientation (no transpose)
        raw_matrix_base = "matrix"

        raw_X = adata.raw.X
        if sp.issparse(raw_X):
            if sp.isspmatrix_csc(raw_X):
                raw_X_to_store = raw_X
            elif sp.isspmatrix_csr(raw_X):
                raw_X_to_store = raw_X.tocsc()
            else:
                raw_X_to_store = sp.csc_matrix(raw_X)
        else:
            raw_X_to_store = sp.csc_matrix(raw_X)

        save_sparse_parquet(raw_X_to_store, raw_dir / raw_matrix_base, orientation='cells_x_genes')
        saved_files['raw_X'] = f"raw/{raw_matrix_base}"
        if compute_hashes:
            hashes['raw_X'] = compute_hash(adata.raw.X)

        # Save raw var (gene metadata for raw counts)
        if adata.raw.var.shape[1] > 0:
            save_parquet(adata.raw.var, raw_dir / "var.parquet")
            saved_files['raw_var'] = "raw/var.parquet"
            if compute_hashes:
                hashes['raw_var'] = compute_hash(adata.raw.var)

        # Save raw gene IDs (always gzip - small file, fast compression)
        raw_features_file = "features.tsv.gz"
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
            compression='gzip'
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
            if compute_hashes:
                hashes['uns'] = compute_hash(uns_json_serializable)

    # =========================================================================
    # Create manifest file
    # =========================================================================
    manifest = {
        'format': 'scio v0.1.2',  # Binary CSC format with cells×genes orientation
        'orientation': 'cells_x_genes',  # v0.1.2: no transpose on Python read, R transposes on read
        'created_by': 'scio (Python)',
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

    # Add hashes if computed
    if compute_hashes and hashes:
        manifest['hashes'] = hashes

    save_json(manifest, folder_path / "manifest.json")

    return saved_files


def update_folder(adata: ad.AnnData,
                  folder_path: Path,
                  manifest: Dict,
                  compress: bool = True) -> Dict[str, int]:
    """
    Update folder with only changed components (incremental update)

    Compares current components with stored hashes and only rewrites
    components that have changed.

    Parameters:
    -----------
    adata : AnnData
        AnnData object to save
    folder_path : Path
        Existing .scio folder path
    manifest : dict
        Existing manifest with hashes
    compress : bool
        Whether to compress files (default: True)

    Returns:
    --------
    dict : {'updated': int, 'skipped': int}
    """
    folder_path = Path(folder_path)
    old_hashes = manifest.get('hashes', {})
    new_hashes = {}
    saved_files = manifest.get('files', {})
    updated_count = 0
    skipped_count = 0

    def check_and_update(component, key: str, save_func):
        nonlocal updated_count, skipped_count
        new_hash = compute_hash(component)
        new_hashes[key] = new_hash

        if old_hashes.get(key) != new_hash:
            save_func()
            updated_count += 1
            print(f"    Updated: {key}", flush=True)
            return True
        else:
            skipped_count += 1
            return False

    # =========================================================================
    # 1. X - Expression matrix (cells×genes, no transpose)
    # =========================================================================
    matrix_base = "matrix"
    saved_files['X'] = matrix_base
    def save_x():
        # v0.1.2: Store in cells×genes orientation (no transpose)
        if sp.issparse(adata.X):
            if sp.isspmatrix_csc(adata.X):
                X_to_store = adata.X
            elif sp.isspmatrix_csr(adata.X):
                X_to_store = adata.X.tocsc()
            else:
                X_to_store = sp.csc_matrix(adata.X)
        else:
            X_to_store = sp.csc_matrix(adata.X)
        save_sparse_parquet(X_to_store, folder_path / matrix_base, orientation='cells_x_genes')
    check_and_update(adata.X, 'X', save_x)

    # =========================================================================
    # 2. obs - Cell metadata
    # =========================================================================
    def save_obs():
        save_parquet(adata.obs, folder_path / "obs.parquet")
        # Also update barcodes (always gzip - small file)
        barcodes_file = "barcodes.tsv.gz"
        barcodes_df = pd.DataFrame(index=adata.obs_names)
        save_tsv_gz(barcodes_df, folder_path / barcodes_file, header=False)
    check_and_update(adata.obs, 'obs', save_obs)

    # =========================================================================
    # 3. var - Gene metadata
    # =========================================================================
    if adata.var.shape[1] > 0:
        def save_var():
            save_parquet(adata.var, folder_path / "var.parquet")
        check_and_update(adata.var, 'var', save_var)

    # =========================================================================
    # 4. obsm - Cell embeddings
    # =========================================================================
    if len(adata.obsm.keys()) > 0:
        obsm_dir = folder_path / "obsm"
        obsm_dir.mkdir(exist_ok=True)

        for key in adata.obsm.keys():
            hash_key = f'obsm_{key}'
            file_path = f"obsm/{key}.parquet"
            saved_files[hash_key] = file_path

            def save_obsm(k=key):
                emb_array = adata.obsm[k]
                emb_df = pd.DataFrame(
                    emb_array,
                    index=adata.obs_names,
                    columns=[f"{k}_{i+1}" for i in range(emb_array.shape[1])]
                )
                save_parquet(emb_df, obsm_dir / f"{k}.parquet")
            check_and_update(adata.obsm[key], hash_key, save_obsm)

    # =========================================================================
    # 5. varm - Gene embeddings
    # =========================================================================
    if len(adata.varm.keys()) > 0:
        varm_dir = folder_path / "varm"
        varm_dir.mkdir(exist_ok=True)

        for key in adata.varm.keys():
            hash_key = f'varm_{key}'
            file_path = f"varm/{key}.parquet"
            saved_files[hash_key] = file_path

            def save_varm(k=key):
                varm_array = adata.varm[k]
                varm_df = pd.DataFrame(
                    varm_array,
                    index=adata.var_names,
                    columns=[f"{k}_{i+1}" for i in range(varm_array.shape[1])]
                )
                save_parquet(varm_df, varm_dir / f"{k}.parquet")
            check_and_update(adata.varm[key], hash_key, save_varm)

    # =========================================================================
    # 6. obsp - Cell-cell graphs
    # =========================================================================
    if len(adata.obsp.keys()) > 0:
        obsp_dir = folder_path / "obsp"
        obsp_dir.mkdir(exist_ok=True)

        for key in adata.obsp.keys():
            hash_key = f'obsp_{key}'
            graph_base = key
            saved_files[hash_key] = f"obsp/{graph_base}"

            def save_obsp(k=key, gb=graph_base):
                save_sparse_parquet(adata.obsp[k], obsp_dir / gb)
            check_and_update(adata.obsp[key], hash_key, save_obsp)

    # =========================================================================
    # 7. varp - Gene-gene graphs
    # =========================================================================
    if len(adata.varp.keys()) > 0:
        varp_dir = folder_path / "varp"
        varp_dir.mkdir(exist_ok=True)

        for key in adata.varp.keys():
            hash_key = f'varp_{key}'
            graph_base = key
            saved_files[hash_key] = f"varp/{graph_base}"

            def save_varp(k=key, gb=graph_base):
                save_sparse_parquet(adata.varp[k], varp_dir / gb)
            check_and_update(adata.varp[key], hash_key, save_varp)

    # =========================================================================
    # 8. layers - Additional matrices
    # =========================================================================
    if len(adata.layers.keys()) > 0:
        layers_dir = folder_path / "layers"
        layers_dir.mkdir(exist_ok=True)

        for key in adata.layers.keys():
            hash_key = f'layer_{key}'
            layer_base = key
            saved_files[hash_key] = f"layers/{layer_base}"

            def save_layer(k=key, lb=layer_base):
                # v0.1.2: Store in cells×genes orientation (no transpose)
                layer = adata.layers[k]
                if sp.issparse(layer):
                    if sp.isspmatrix_csc(layer):
                        layer_to_store = layer
                    elif sp.isspmatrix_csr(layer):
                        layer_to_store = layer.tocsc()
                    else:
                        layer_to_store = sp.csc_matrix(layer)
                else:
                    layer_to_store = sp.csc_matrix(layer)
                save_sparse_parquet(layer_to_store, layers_dir / lb, orientation='cells_x_genes')
            check_and_update(adata.layers[key], hash_key, save_layer)

    # =========================================================================
    # 9. raw - Raw data
    # =========================================================================
    if adata.raw is not None:
        raw_dir = folder_path / "raw"
        raw_dir.mkdir(exist_ok=True)

        raw_matrix_base = "matrix"
        saved_files['raw_X'] = f"raw/{raw_matrix_base}"
        def save_raw_x():
            # v0.1.2: Store in cells×genes orientation (no transpose)
            raw_X = adata.raw.X
            if sp.issparse(raw_X):
                if sp.isspmatrix_csc(raw_X):
                    raw_X_to_store = raw_X
                elif sp.isspmatrix_csr(raw_X):
                    raw_X_to_store = raw_X.tocsc()
                else:
                    raw_X_to_store = sp.csc_matrix(raw_X)
            else:
                raw_X_to_store = sp.csc_matrix(raw_X)
            save_sparse_parquet(raw_X_to_store, raw_dir / raw_matrix_base, orientation='cells_x_genes')
        check_and_update(adata.raw.X, 'raw_X', save_raw_x)

        if adata.raw.var.shape[1] > 0:
            def save_raw_var():
                save_parquet(adata.raw.var, raw_dir / "var.parquet")
            check_and_update(adata.raw.var, 'raw_var', save_raw_var)

    # =========================================================================
    # 10. uns - Unstructured metadata
    # =========================================================================
    if len(adata.uns.keys()) > 0:
        uns_json_serializable = {}
        for key, value in adata.uns.items():
            try:
                serializable_value = numpy_to_json_serializable(value)
                import json
                json.dumps(serializable_value)
                uns_json_serializable[key] = serializable_value
            except (TypeError, ValueError):
                warnings.warn(f"Skipping non-serializable uns entry: {key}")
                continue

        if uns_json_serializable:
            def save_uns():
                save_json(uns_json_serializable, folder_path / "uns.json")
            check_and_update(uns_json_serializable, 'uns', save_uns)

    # =========================================================================
    # Update manifest
    # =========================================================================
    manifest['format'] = 'scio v0.1.2'
    manifest['orientation'] = 'cells_x_genes'
    manifest['hashes'] = new_hashes
    manifest['files'] = saved_files

    # Update components list
    manifest['components'] = {
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
    }

    save_json(manifest, folder_path / "manifest.json")

    print(f"    Summary: {updated_count} components updated, {skipped_count} unchanged", flush=True)

    return {'updated': updated_count, 'skipped': skipped_count}
