
import os
import pandas as pd
import scipy.io
import anndata as ad
import memento

from scipy.sparse import csr_matrix

# -------------------------
# Base directory
# -------------------------
data_path = r"/project/gilad/awchen55/differentialDispersion/data/simulations/mean_var_independence/"

# -------------------------
# File paths
# -------------------------
counts_path   = os.path.join(data_path, "counts.mtx")
genes_path    = os.path.join(data_path, "genes.txt")
barcodes_path = os.path.join(data_path, "barcodes.txt")
metadata_path = os.path.join(data_path, "cell_metadata.csv")

# -------------------------
# Load counts matrix
# Matrix Market is often genes x cells
# AnnData requires cells x genes
# -------------------------
counts = scipy.io.mmread(counts_path).tocsr()

# -------------------------
# Load genes
# Adjust sep if needed
# -------------------------
genes = pd.read_csv(genes_path, header=None, sep="\t")
genes.columns = ["gene"]

# -------------------------
# Load barcodes
# -------------------------
barcodes = pd.read_csv(barcodes_path, header=None)
barcodes.columns = ["barcode"]

# -------------------------
# Load cell metadata
# -------------------------
metadata = pd.read_csv(metadata_path)

# -------------------------
# Basic dimension checks
# -------------------------
print("counts shape (raw):", counts.shape)
print("number of genes:", genes.shape[0])
print("number of barcodes:", barcodes.shape[0])
print("metadata shape:", metadata.shape)

# Check whether counts is genes x cells
if counts.shape[0] == genes.shape[0] and counts.shape[1] == barcodes.shape[0]:
    print("Detected counts as genes x cells; transposing for AnnData.")
    X = counts.T.tocsr()
elif counts.shape[0] == barcodes.shape[0] and counts.shape[1] == genes.shape[0]:
    print("Detected counts as cells x genes; no transpose needed.")
    X = counts.tocsr()
else:
    raise ValueError(
        "Counts dimensions do not match genes/barcodes.\n"
        f"counts shape: {counts.shape}\n"
        f"genes: {genes.shape[0]}\n"
        f"barcodes: {barcodes.shape[0]}"
    )

# -------------------------
# Align metadata to barcodes
# Assumes metadata rows are in same order as barcodes
# -------------------------
if metadata.shape[0] != barcodes.shape[0]:
    raise ValueError(
        "Number of rows in cell_metadata.csv does not match number of barcodes.\n"
        f"metadata rows: {metadata.shape[0]}\n"
        f"barcodes: {barcodes.shape[0]}"
    )

metadata.index = barcodes["barcode"].astype(str).values

# -------------------------
# Create gene metadata
# -------------------------
genes["gene"] = genes["gene"].astype(str)
genes.index = genes["gene"].values

# -------------------------
# Create AnnData object
# -------------------------
adata = ad.AnnData(
    X=X,
    obs=metadata,
    var=genes
)

adata.obs_names = barcodes["barcode"].astype(str).values
adata.var_names = genes["gene"].astype(str).values

# Make names unique just in case
adata.obs_names_make_unique()
adata.var_names_make_unique()

# Force CSR matrix for memento
adata.X = csr_matrix(adata.X)

print(adata)
print("adata.X type:", type(adata.X))

# -------------------------
# Save AnnData object
# -------------------------
h5ad_path = os.path.join(data_path, "sim_cardiomyocytes.h5ad")
adata.write(h5ad_path)
print(f"Saved AnnData to: {h5ad_path}")

# -------------------------
# Setup memento
# -------------------------
adata.obs["capture_rate"] = 0.2

memento.setup_memento(adata, q_column="capture_rate")
memento.create_groups(
    adata,
    label_columns=["individual", "replicate", "celltype"]
)

# -------------------------
# Compute 1D moments
# -------------------------
memento.compute_1d_moments(
    adata,
    min_perc_group=0.95
)

# -------------------------
# Extract mean / variance / counts
# -------------------------
mean, var, counts_out = memento.get_1d_moments(adata)

# -------------------------
# Save outputs
# -------------------------
var_path   = os.path.join(data_path, "memento_var_estimate_simulation.csv")
mean_path  = os.path.join(data_path, "memento_mean_estimate_simulation.csv")
count_path = os.path.join(data_path, "memento_counts_simulation.csv")

var.to_csv(var_path)
mean.to_csv(mean_path)
counts_out.to_csv(count_path)

print(f"Saved variance estimates to: {var_path}")
print(f"Saved mean estimates to: {mean_path}")
print(f"Saved counts to: {count_path}")
    




