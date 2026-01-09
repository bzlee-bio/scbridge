"""
Lazy loading module - Load AnnData components on-demand for better performance
"""
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, Any, Dict

from .formats import load_json, load_parquet
from .readers import (
    _detect_gzip,
    _ensure_list,
    _load_sparse_matrix
)


class LazyAnnData:
    """
    Lazy-loading wrapper for AnnData that only loads components when accessed.

    This significantly improves performance for large datasets by:
    - Loading only manifest and basic metadata initially
    - Loading data components (X, obs, obsm, etc.) on first access
    - Caching loaded components in memory for subsequent access

    Example:
    --------
    >>> import scio
    >>> # Fast - only loads structure
    >>> adata = scio.read("data.scio", lazy=True)
    >>>
    >>> # Instant - just metadata
    >>> print(adata.shape)  # (300000, 20000)
    >>>
    >>> # Loads obs on first access
    >>> cell_types = adata.obs['cell_type']
    >>>
    >>> # Loads X on first access
    >>> expression = adata.X
    """

    def __init__(self, folder_path: Path):
        """
        Initialize LazyAnnData from a folder path

        Parameters:
        -----------
        folder_path : Path
            Path to .scio folder
        """
        self._folder_path = Path(folder_path)
        self._cache: Dict[str, Any] = {}

        # Load manifest immediately (small file)
        manifest_path = self._folder_path / "manifest.json"
        if not manifest_path.exists():
            raise FileNotFoundError(f"Manifest file not found: {manifest_path}")

        self._manifest = load_json(manifest_path)

        # Load obs_names and var_names immediately (small, needed for shape)
        self._load_names()

        # Mark what's not loaded yet
        self._loaded = {
            'X': False,
            'obs': False,
            'var': False,
            'obsm': False,
            'varm': False,
            'obsp': False,
            'varp': False,
            'layers': False,
            'uns': False,
            'raw': False
        }

    def _load_names(self):
        """Load obs_names and var_names (small files, needed for shape)"""
        # Load barcodes (obs_names)
        barcodes_file = self._folder_path / self._manifest['files']['barcodes']
        compression = 'gzip' if _detect_gzip(barcodes_file) else None
        barcodes_df = pd.read_csv(barcodes_file, sep='\t', header=None, compression=compression)
        self._cache['obs_names'] = barcodes_df.iloc[:, 0].values

        # Load features (var_names)
        features_file = self._folder_path / self._manifest['files']['features']
        compression = 'gzip' if _detect_gzip(features_file) else None
        features_df = pd.read_csv(features_file, sep='\t', header=None, compression=compression)
        self._cache['var_names'] = features_df.iloc[:, 0].values

    def _load_X(self):
        """Lazy load X matrix"""
        if not self._loaded['X']:
            self._cache['X'] = _load_sparse_matrix(
                self._folder_path, 'X', self._manifest, is_expression_matrix=True
            )
            self._loaded['X'] = True
        return self._cache['X']

    def _load_obs(self):
        """Lazy load obs DataFrame"""
        if not self._loaded['obs']:
            obs_file = self._folder_path / self._manifest['files']['obs']
            self._cache['obs'] = load_parquet(obs_file)
            self._loaded['obs'] = True
        return self._cache['obs']

    def _load_var(self):
        """Lazy load var DataFrame"""
        if not self._loaded['var']:
            if 'var' in self._manifest['files']:
                var_file = self._folder_path / self._manifest['files']['var']
                self._cache['var'] = load_parquet(var_file)
            else:
                # Create empty var with correct index
                self._cache['var'] = pd.DataFrame(index=self._cache['var_names'])
            self._loaded['var'] = True
        return self._cache['var']

    def _load_obsm(self):
        """Lazy load obsm (cell embeddings)"""
        if not self._loaded['obsm']:
            obsm = {}
            obsm_keys = _ensure_list(self._manifest['components']['obsm'])
            for key in obsm_keys:
                file_path = self._folder_path / self._manifest['files'][f'obsm_{key}']
                emb_df = load_parquet(file_path)
                obsm[key] = emb_df.values
            self._cache['obsm'] = obsm
            self._loaded['obsm'] = True
        return self._cache['obsm']

    def _load_varm(self):
        """Lazy load varm (gene embeddings)"""
        if not self._loaded['varm']:
            varm = {}
            varm_keys = _ensure_list(self._manifest['components']['varm'])
            for key in varm_keys:
                file_path = self._folder_path / self._manifest['files'][f'varm_{key}']
                varm_df = load_parquet(file_path)
                varm[key] = varm_df.values
            self._cache['varm'] = varm
            self._loaded['varm'] = True
        return self._cache['varm']

    def _load_obsp(self):
        """Lazy load obsp (cell-cell graphs)"""
        if not self._loaded['obsp']:
            obsp = {}
            obsp_keys = _ensure_list(self._manifest['components']['obsp'])
            for key in obsp_keys:
                obsp[key] = _load_sparse_matrix(
                    self._folder_path, f'obsp_{key}', self._manifest, transpose=False
                )
            self._cache['obsp'] = obsp
            self._loaded['obsp'] = True
        return self._cache['obsp']

    def _load_varp(self):
        """Lazy load varp (gene-gene graphs)"""
        if not self._loaded['varp']:
            varp = {}
            varp_keys = _ensure_list(self._manifest['components']['varp'])
            for key in varp_keys:
                varp[key] = _load_sparse_matrix(
                    self._folder_path, f'varp_{key}', self._manifest, transpose=False
                )
            self._cache['varp'] = varp
            self._loaded['varp'] = True
        return self._cache['varp']

    def _load_layers(self):
        """Lazy load layers (additional matrices)"""
        if not self._loaded['layers']:
            layers = {}
            layer_keys = _ensure_list(self._manifest['components']['layers'])
            for key in layer_keys:
                layers[key] = _load_sparse_matrix(
                    self._folder_path, f'layer_{key}', self._manifest, is_expression_matrix=True
                )
            self._cache['layers'] = layers
            self._loaded['layers'] = True
        return self._cache['layers']

    def _load_uns(self):
        """Lazy load uns (unstructured metadata)"""
        if not self._loaded['uns']:
            if 'uns' in self._manifest['files']:
                uns_file = self._folder_path / self._manifest['files']['uns']
                self._cache['uns'] = load_json(uns_file)
            else:
                self._cache['uns'] = {}
            self._loaded['uns'] = True
        return self._cache['uns']

    def _load_raw(self):
        """Lazy load raw data"""
        if not self._loaded['raw']:
            if self._manifest['components']['raw']:
                raw_dir = self._folder_path / "raw"

                # Load raw X matrix
                raw_X = _load_sparse_matrix(
                    self._folder_path, 'raw_X', self._manifest, is_expression_matrix=True
                )

                # Load raw gene IDs
                raw_features_file = raw_dir / self._manifest['files']['raw_features'].split('/')[-1]
                compression = 'gzip' if _detect_gzip(raw_features_file) else None
                raw_features_df = pd.read_csv(
                    raw_features_file,
                    sep='\t',
                    header=None,
                    compression=compression
                )
                raw_var_names = raw_features_df.iloc[:, 0].values

                # Load raw var metadata if present
                if 'raw_var' in self._manifest['files']:
                    raw_var_file = self._folder_path / self._manifest['files']['raw_var']
                    raw_var = load_parquet(raw_var_file)
                else:
                    raw_var = pd.DataFrame(index=raw_var_names)

                # Create raw AnnData
                raw_adata = ad.AnnData(X=raw_X, var=raw_var)
                raw_adata.var_names = raw_var_names
                raw_adata.obs_names = self._cache['obs_names']

                self._cache['raw'] = raw_adata
            else:
                self._cache['raw'] = None
            self._loaded['raw'] = True
        return self._cache['raw']

    # Property accessors for AnnData-like interface
    @property
    def X(self):
        """Access X matrix (loads on first access)"""
        return self._load_X()

    @property
    def obs(self):
        """Access obs DataFrame (loads on first access)"""
        return self._load_obs()

    @property
    def var(self):
        """Access var DataFrame (loads on first access)"""
        return self._load_var()

    @property
    def obsm(self):
        """Access obsm dict (loads on first access)"""
        return self._load_obsm()

    @property
    def varm(self):
        """Access varm dict (loads on first access)"""
        return self._load_varm()

    @property
    def obsp(self):
        """Access obsp dict (loads on first access)"""
        return self._load_obsp()

    @property
    def varp(self):
        """Access varp dict (loads on first access)"""
        return self._load_varp()

    @property
    def layers(self):
        """Access layers dict (loads on first access)"""
        return self._load_layers()

    @property
    def uns(self):
        """Access uns dict (loads on first access)"""
        return self._load_uns()

    @property
    def raw(self):
        """Access raw data (loads on first access)"""
        return self._load_raw()

    @property
    def obs_names(self):
        """Get obs_names (already loaded)"""
        return self._cache['obs_names']

    @property
    def var_names(self):
        """Get var_names (already loaded)"""
        return self._cache['var_names']

    @property
    def shape(self):
        """Get shape without loading data"""
        return (len(self._cache['obs_names']), len(self._cache['var_names']))

    @property
    def n_obs(self):
        """Number of observations"""
        return len(self._cache['obs_names'])

    @property
    def n_vars(self):
        """Number of variables"""
        return len(self._cache['var_names'])

    def to_memory(self) -> ad.AnnData:
        """
        Load all components into memory and return a regular AnnData object.

        This is useful when you want to convert from lazy loading to full loading,
        for example before performing operations that require the full dataset.

        Returns:
        --------
        AnnData : Fully loaded AnnData object
        """
        # Load all components
        X = self._load_X()
        obs = self._load_obs()
        var = self._load_var()

        # Create base AnnData
        adata = ad.AnnData(X=X, obs=obs, var=var)
        adata.obs_names = self._cache['obs_names']
        adata.var_names = self._cache['var_names']

        # Load optional components
        obsm = self._load_obsm()
        for key, val in obsm.items():
            adata.obsm[key] = val

        varm = self._load_varm()
        for key, val in varm.items():
            adata.varm[key] = val

        obsp = self._load_obsp()
        for key, val in obsp.items():
            adata.obsp[key] = val

        varp = self._load_varp()
        for key, val in varp.items():
            adata.varp[key] = val

        layers = self._load_layers()
        for key, val in layers.items():
            adata.layers[key] = val

        adata.uns = self._load_uns()
        adata.raw = self._load_raw()

        return adata

    def __repr__(self):
        """String representation"""
        loaded_components = [k for k, v in self._loaded.items() if v]
        loaded_str = f" (loaded: {', '.join(loaded_components)})" if loaded_components else " (no components loaded yet)"
        return f"LazyAnnData object with n_obs × n_vars = {self.n_obs} × {self.n_vars}{loaded_str}"

    def __str__(self):
        """String representation"""
        return self.__repr__()
