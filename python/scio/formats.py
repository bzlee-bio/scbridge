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
        # For sparse matrices, hash the data, indices, and shape
        data = (obj.data.tobytes(), obj.indices.tobytes(),
                obj.indptr.tobytes() if hasattr(obj, 'indptr') else b'',
                str(obj.shape).encode())
        return hashlib.md5(b''.join(data)).hexdigest()
    elif isinstance(obj, np.ndarray):
        return hashlib.md5(obj.tobytes()).hexdigest()
    elif isinstance(obj, pd.DataFrame):
        # Hash DataFrame contents using a stable representation
        # Convert to parquet bytes for consistent hashing
        import io
        buffer = io.BytesIO()
        obj.to_parquet(buffer)
        return hashlib.md5(buffer.getvalue()).hexdigest()
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
