#!/usr/bin/env python3
"""
Test script for scio Python package
Tests write/read round-trip with AnnData objects
"""

import sys
import tempfile
import shutil
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
    import anndata as ad
    from scipy import sparse
    import scio
except ImportError as e:
    print(f"Error importing required packages: {e}")
    print("Please install: pip install numpy pandas anndata scipy")
    sys.exit(1)

def create_test_adata():
    """Create a simple test AnnData object"""
    np.random.seed(42)

    # Create count matrix (20 genes x 10 cells)
    counts = sparse.random(20, 10, density=0.3, format='csr')
    counts.data = np.round(counts.data * 10)

    # Create AnnData
    adata = ad.AnnData(X=counts)

    # Add obs (cell metadata)
    adata.obs['sample'] = pd.Categorical(['A'] * 5 + ['B'] * 5)
    adata.obs['celltype'] = pd.Categorical(['T', 'B', 'NK'] * 3 + ['T'])

    # Add var (gene metadata)
    adata.var['gene_name'] = [f'Gene{i}' for i in range(20)]
    adata.var['highly_variable'] = np.random.choice([True, False], 20)

    # Add obsm (embeddings)
    adata.obsm['X_pca'] = np.random.randn(10, 5)
    adata.obsm['X_umap'] = np.random.randn(10, 2)

    # Add layers
    adata.layers['normalized'] = counts.copy()
    adata.layers['scaled'] = counts.copy()

    # Add uns (unstructured)
    adata.uns['method'] = 'test'
    adata.uns['params'] = {'n_neighbors': 15}

    return adata

def test_write_read():
    """Test write and read functionality"""
    print("=== Creating test AnnData ===")
    adata = create_test_adata()
    print(f"Original AnnData: {adata}")
    print(f"  Shape: {adata.shape}")
    print(f"  Layers: {list(adata.layers.keys())}")
    print(f"  Obsm: {list(adata.obsm.keys())}")

    # Create temporary directory for test
    test_dir = Path(tempfile.mkdtemp())
    test_path = test_dir / "test_data.scio"

    try:
        print(f"\n=== Testing scio.write() ===")
        print(f"Writing to: {test_path}")
        scio.write(adata, str(test_path))

        # Check what files were created
        print("\n=== Checking created files ===")
        if test_path.exists():
            files = list(test_path.rglob("*"))
            files = [f.relative_to(test_path) for f in files if f.is_file()]
            print("Files created:")
            for f in sorted(files):
                size = (test_path / f).stat().st_size / 1024
                print(f"  {f} ({size:.2f} KB)")
        else:
            print("ERROR: Output directory not created!")
            return False

        print(f"\n=== Testing scio.read() ===")
        print(f"Reading from: {test_path}")
        adata_back = scio.read(str(test_path))

        print(f"\nLoaded AnnData: {adata_back}")
        print(f"  Shape: {adata_back.shape}")
        print(f"  Layers: {list(adata_back.layers.keys())}")
        print(f"  Obsm: {list(adata_back.obsm.keys())}")

        # Verify data integrity
        print("\n=== Verifying data integrity ===")
        checks = {
            "Shape matches": adata.shape == adata_back.shape,
            "X matrix matches": np.allclose(adata.X.toarray(), adata_back.X.toarray()),
            "obs matches": adata.obs.equals(adata_back.obs),
            "var shape matches": adata.var.shape == adata_back.var.shape,
            "Layers match": list(adata.layers.keys()) == list(adata_back.layers.keys()),
            "Obsm keys match": list(adata.obsm.keys()) == list(adata_back.obsm.keys()),
            "PCA embeddings match": np.allclose(adata.obsm['X_pca'], adata_back.obsm['X_pca']),
            "Uns preserved": 'method' in adata_back.uns
        }

        for check_name, result in checks.items():
            status = "✓" if result else "✗"
            print(f"{status} {check_name}")

        all_passed = all(checks.values())

        if all_passed:
            print("\n=== ✓ ALL PYTHON TESTS PASSED ===")
        else:
            print("\n=== ✗ SOME PYTHON TESTS FAILED ===")

        return all_passed

    finally:
        # Clean up
        print("\n=== Cleaning up ===")
        if test_dir.exists():
            shutil.rmtree(test_dir)
            print("Test directory removed")

if __name__ == "__main__":
    success = test_write_read()
    sys.exit(0 if success else 1)
