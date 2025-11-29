"""
Format handlers for different file types (MTX, Parquet, JSON)
"""
import pandas as pd
import numpy as np
import scipy.sparse as sp
import scipy.io as sio
import json
import gzip
import hashlib
import pickle
from pathlib import Path
from typing import Union, Dict, Any


def compute_hash(obj: Any) -> str:
    """
    Compute MD5 hash of a Python object for change detection.

    Used for incremental updates - only rewrite components that changed.

    IMPORTANT: This function produces STABLE hashes that are independent of:
    - Memory layout (C vs Fortran order for numpy arrays)
    - Index type (RangeIndex vs Int64Index for DataFrames)
    - Internal pandas/numpy representation details

    Parameters:
    -----------
    obj : Any
        Object to hash (numpy array, sparse matrix, DataFrame, dict, etc.)

    Returns:
    --------
    str : MD5 hash string
    """
    # Convert object to bytes for hashing
    if sp.issparse(obj):
        # For sparse matrices, hash internal CSR/CSC arrays directly
        # This is much faster than converting to COO for large matrices
        if sp.isspmatrix_csr(obj):
            # CSR: hash data, indices, indptr directly
            hasher = hashlib.md5()
            hasher.update(obj.data.tobytes())
            hasher.update(obj.indices.tobytes())
            hasher.update(obj.indptr.tobytes())
            hasher.update(str(obj.shape).encode())
            return hasher.hexdigest()
        elif sp.isspmatrix_csc(obj):
            # CSC: hash data, indices, indptr directly
            hasher = hashlib.md5()
            hasher.update(obj.data.tobytes())
            hasher.update(obj.indices.tobytes())
            hasher.update(obj.indptr.tobytes())
            hasher.update(str(obj.shape).encode())
            return hasher.hexdigest()
        else:
            # For other sparse formats, convert to CSR first
            csr = obj.tocsr()
            hasher = hashlib.md5()
            hasher.update(csr.data.tobytes())
            hasher.update(csr.indices.tobytes())
            hasher.update(csr.indptr.tobytes())
            hasher.update(str(csr.shape).encode())
            return hasher.hexdigest()
    elif isinstance(obj, np.ndarray):
        # Use contiguous array to ensure consistent byte representation
        # regardless of memory layout (C vs Fortran order)
        arr = np.ascontiguousarray(obj)
        return hashlib.md5(arr.tobytes()).hexdigest()
    elif isinstance(obj, pd.DataFrame):
        # Hash DataFrame contents using a stable representation
        # Reset index to avoid RangeIndex vs Int64Index differences
        # Sort columns for consistency
        df = obj.reset_index(drop=True)
        # Use values + column names + dtypes for stable hash
        hasher = hashlib.md5()
        hasher.update(str(sorted(df.columns.tolist())).encode())
        hasher.update(str(df.dtypes.tolist()).encode())
        for col in sorted(df.columns):
            col_data = df[col]
            if col_data.dtype == 'category':
                # For categorical, hash codes and categories
                hasher.update(col_data.cat.codes.values.tobytes())
                hasher.update(str(col_data.cat.categories.tolist()).encode())
            elif col_data.dtype == 'object':
                # For object dtype (strings), convert to string and encode
                # This avoids memory address differences
                hasher.update('|'.join(col_data.astype(str).values).encode('utf-8'))
            elif col_data.dtype == 'bool':
                # Bool arrays need explicit conversion
                hasher.update(col_data.values.astype(np.uint8).tobytes())
            else:
                hasher.update(np.ascontiguousarray(col_data.values).tobytes())
        return hasher.hexdigest()
    else:
        # For other objects, use pickle
        return hashlib.md5(pickle.dumps(obj)).hexdigest()


def save_mtx(matrix: Union[np.ndarray, sp.spmatrix],
             path: Path,
             compress: bool = True) -> Path:
    """
    Save a matrix in Matrix Market (MTX) format

    Parameters:
    -----------
    matrix : array-like
        Matrix to save (will be converted to sparse if needed)
    path : Path
        Output file path
    compress : bool
        Whether to gzip compress (default: True)

    Returns:
    --------
    Path : Path to saved file
    """
    # Convert to sparse if needed
    if not sp.issparse(matrix):
        matrix = sp.csr_matrix(matrix)

    # Add .gz extension if compressing
    if compress and not str(path).endswith('.gz'):
        path = Path(str(path) + '.gz')

    # Save
    if compress:
        with gzip.open(path, 'wb') as f:
            sio.mmwrite(f, matrix)
    else:
        sio.mmwrite(path, matrix)

    return path


def load_mtx(path: Path, transpose: bool = False) -> sp.csr_matrix:
    """
    Load a matrix from Matrix Market (MTX) format

    Parameters:
    -----------
    path : Path
        Path to MTX file (can be gzipped)
    transpose : bool
        Whether to transpose after loading (default: False)

    Returns:
    --------
    csr_matrix : Loaded sparse matrix
    """
    # Auto-detect gzip
    if str(path).endswith('.gz'):
        with gzip.open(path, 'rb') as f:
            matrix = sio.mmread(f)
    else:
        matrix = sio.mmread(path)

    # Convert to CSR format
    matrix = sp.csr_matrix(matrix)

    if transpose:
        matrix = matrix.T

    return matrix


def save_parquet(df: pd.DataFrame, path: Path) -> Path:
    """
    Save a DataFrame to Parquet format

    Parameters:
    -----------
    df : DataFrame
        DataFrame to save
    path : Path
        Output file path

    Returns:
    --------
    Path : Path to saved file
    """
    df.to_parquet(path)
    return path


def load_parquet(path: Path) -> pd.DataFrame:
    """
    Load a DataFrame from Parquet format

    Parameters:
    -----------
    path : Path
        Path to Parquet file

    Returns:
    --------
    DataFrame : Loaded DataFrame
    """
    return pd.read_parquet(path)


def save_json(data: Dict[str, Any], path: Path, indent: int = 2) -> Path:
    """
    Save a dictionary to JSON format

    Parameters:
    -----------
    data : dict
        Dictionary to save
    path : Path
        Output file path
    indent : int
        JSON indentation level (default: 2)

    Returns:
    --------
    Path : Path to saved file
    """
    with open(path, 'w') as f:
        json.dump(data, f, indent=indent)
    return path


def load_json(path: Path) -> Dict[str, Any]:
    """
    Load a dictionary from JSON format

    Parameters:
    -----------
    path : Path
        Path to JSON file

    Returns:
    --------
    dict : Loaded dictionary
    """
    with open(path, 'r') as f:
        return json.load(f)


def save_tsv_gz(data: pd.DataFrame, path: Path, header: bool = False) -> Path:
    """
    Save data to gzipped TSV format

    Parameters:
    -----------
    data : DataFrame
        Data to save
    path : Path
        Output file path
    header : bool
        Whether to include header (default: False)

    Returns:
    --------
    Path : Path to saved file
    """
    data.to_csv(path, sep='\t', header=header, compression='gzip')
    return path


def load_tsv_gz(path: Path, header: bool = False) -> pd.DataFrame:
    """
    Load data from gzipped TSV format

    Parameters:
    -----------
    path : Path
        Path to TSV file
    header : bool
        Whether file has header (default: False)

    Returns:
    --------
    DataFrame : Loaded data
    """
    if header:
        return pd.read_csv(path, sep='\t', compression='gzip')
    else:
        return pd.read_csv(path, sep='\t', compression='gzip', header=None)


def save_sparse_parquet(matrix: Union[np.ndarray, sp.spmatrix],
                        path: Path,
                        shape_file: Path = None,
                        orientation: str = None) -> Path:
    """
    Save a sparse matrix in binary format (CSC or CSR).

    Supports both CSC and CSR formats to avoid expensive conversion.
    This is the fastest format - directly saves sparse components as numpy binary files.
    Much faster than MTX and Parquet for large sparse matrices.

    Files created:
    - {path}.data.npy - non-zero values
    - {path}.indices.npy - row/column indices
    - {path}.indptr.npy - column/row pointers
    - {path}.shape.json - metadata (shape, dtype, nnz, format, orientation)

    Parameters:
    -----------
    matrix : array-like
        Matrix to save (CSC or CSR - saved directly without conversion)
    path : Path
        Base output file path (without extension)
    shape_file : Path, optional
        Path for shape metadata JSON file. If None, uses {path}.shape.json
    orientation : str, optional
        Data orientation: 'cells_x_genes' (v0.1.2+) or 'genes_x_cells' (legacy).
        Stored in metadata to help readers decide if transpose is needed.

    Returns:
    --------
    Path : Base path to saved files
    """
    path = Path(path)

    # Support both CSC and CSR formats to avoid expensive conversion
    if not sp.issparse(matrix):
        matrix = sp.csc_matrix(matrix)
        sparse_format = 'csc_binary'
    elif sp.isspmatrix_csc(matrix):
        sparse_format = 'csc_binary'
    elif sp.isspmatrix_csr(matrix):
        sparse_format = 'csr_binary'
    else:
        # Convert other formats (COO, etc.) to CSC
        matrix = matrix.tocsc()
        sparse_format = 'csc_binary'

    # Save sparse components as binary numpy files in parallel for faster I/O
    import concurrent.futures

    def _save_npy(filepath, arr):
        np.save(filepath, arr)

    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        futures = [
            executor.submit(_save_npy, str(path) + '.data.npy', matrix.data),
            executor.submit(_save_npy, str(path) + '.indices.npy', matrix.indices),
            executor.submit(_save_npy, str(path) + '.indptr.npy', matrix.indptr),
        ]
        # Wait for all to complete
        for f in futures:
            f.result()

    # Save shape metadata
    if shape_file is None:
        shape_file = Path(str(path) + '.shape.json')

    shape_meta = {
        'shape': list(matrix.shape),
        'dtype': str(matrix.dtype),
        'nnz': matrix.nnz,
        'format': sparse_format
    }

    # Add orientation metadata for v0.1.2+
    if orientation:
        shape_meta['orientation'] = orientation

    with open(shape_file, 'w') as f:
        json.dump(shape_meta, f)

    return path


def load_sparse_parquet(path: Path,
                        shape_file: Path = None,
                        transpose: bool = False) -> sp.spmatrix:
    """
    Load a sparse matrix from binary CSC format.

    When transpose=True, returns CSR format directly (zero-copy from CSC.T).
    This is the key optimization: CSC.T gives CSR as a view, no conversion needed.

    Parameters:
    -----------
    path : Path
        Base path to binary files (without .data.npy extension)
    shape_file : Path, optional
        Path to shape metadata JSON file. If None, uses {path}.shape.json
    transpose : bool
        Whether to transpose after loading (default: False).
        If True, returns CSR matrix. If False, returns CSC matrix.

    Returns:
    --------
    spmatrix : Loaded sparse matrix (CSR if transposed, CSC otherwise)
    """
    path = Path(path)

    # Load shape metadata
    if shape_file is None:
        shape_file = Path(str(path) + '.shape.json')

    with open(shape_file, 'r') as f:
        shape_meta = json.load(f)

    shape = tuple(shape_meta['shape'])
    dtype = np.dtype(shape_meta['dtype'])
    stored_format = shape_meta.get('format', 'csc_binary')

    # Load components from binary files in parallel for faster I/O
    import concurrent.futures

    def _load_npy(filepath):
        return np.load(filepath)

    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        data_future = executor.submit(_load_npy, str(path) + '.data.npy')
        indices_future = executor.submit(_load_npy, str(path) + '.indices.npy')
        indptr_future = executor.submit(_load_npy, str(path) + '.indptr.npy')

        data = data_future.result()
        indices = indices_future.result()
        indptr = indptr_future.result()

    # Reconstruct matrix based on stored format
    # Only convert dtype if necessary (avoid expensive copy for same dtype)
    if data.dtype != dtype:
        data = data.astype(dtype)

    if stored_format == 'csc_binary':
        # New v2.0 CSC format - use copy=False to avoid expensive validation
        matrix = sp.csc_matrix((data, indices, indptr), shape=shape, copy=False)
        if transpose:
            # CSC.T gives CSR directly (zero-copy view) - this is the key optimization!
            matrix = matrix.T
    else:
        # Legacy CSR format (for backward compatibility)
        matrix = sp.csr_matrix((data, indices, indptr), shape=shape, copy=False)
        if transpose:
            # CSR.T gives CSC, need to convert back to CSR (slow)
            matrix = matrix.T.tocsr()

    return matrix


def numpy_to_json_serializable(obj: Any) -> Any:
    """
    Convert numpy types to JSON-serializable types

    Parameters:
    -----------
    obj : Any
        Object to convert

    Returns:
    --------
    Any : JSON-serializable version
    """
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.bool_):
        return bool(obj)
    elif isinstance(obj, dict):
        return {k: numpy_to_json_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [numpy_to_json_serializable(item) for item in obj]
    else:
        return obj
