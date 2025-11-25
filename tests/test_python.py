import scio
import anndata
import scanpy as sc
adata = sc.read_h5ad("/mnt/dks_nas2/sce_qc.h5ad")
# adata = sc.read_h5ad("/mnt/gluster_server/adata_processed.combined.h5ad")
scio.write(adata, "./sce_qc.scio")