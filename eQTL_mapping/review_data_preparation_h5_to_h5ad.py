import h5py
import numpy as np
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import pandas as pd

for treatment in ["IFN", "LPS", "untreated"]:

    file_path = f"../../data/for_saigeQTL/{treatment}.h5"
    # Inspect HDF5 structure
    with h5py.File(file_path, "r") as f:
        def print_structure(name, obj):
            print(name, type(obj))
        f.visititems(print_structure)
    # Read AnnData
    adata = sc.read_10x_h5(file_path)
    # Load metadata
    meta = pd.read_csv(
        f"../../data/for_saigeQTL/{treatment}_metadata.csv",
        index_col=0
    )
    # Check that cell numbers match
    print(meta.shape)
    print(adata.n_obs)
    # Reorder metadata to match adata
    meta = meta.loc[adata.obs_names]
    # Add metadata
    adata.obs = meta
    # Save
    adata.write_h5ad(
        f"../../data/for_saigeQTL/{treatment}.h5ad"
    )

    
### save subsets for comparison with tensorQTL

import h5py
import numpy as np
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import pandas as pd

adata = sc.read_h5ad("../../data/for_saigeQTL/untreated.h5ad")


meta_tensorqtl = pd.read_csv(
    "../../data/for_tensorQTL/tensorQTL_metadata_sum_sizefactorsNorm_log2_untreated_Not_proliferating.txt",
    sep="\t"
)
keep_lines = meta_tensorqtl.columns.tolist()
keep_lines = keep_lines[0].split()
keep_lines = [x.replace(".", "-") for x in keep_lines]

adata_sub = adata[
    (adata.obs["donor_id"].isin(keep_lines)),
].copy()

print("cells:", adata_sub.n_obs)

## dummy because otherwise QTLight throws error
adata_sub.obs["Azimuth:predicted.celltype.l1"] = "all"
print(adata_sub.obs.head())

genes = pd.read_csv(
    "../../../resources/biomart/Homo_sapiens.GRCh38.111.genes.csv"
)

print(genes.head())
name_to_id = genes.drop_duplicates("gene_name").set_index("gene_name")["gene_id"]
mapped = adata_sub.var_names.isin(name_to_id.index)

adata_sub.var["gene_name"] = adata_sub.var_names

adata_sub.var["gene_id"] = (
    adata_sub.var_names
    .map(name_to_id)
)

# Keep only genes that mapped
adata_sub = adata_sub[:, ~adata_sub.var["gene_id"].isna()].copy()

# Make gene_id the primary identifier
adata_sub.var_names = adata_sub.var["gene_id"]

print(f"{mapped.sum()} / {adata_sub.n_vars} genes mapped")
print(adata_sub.var.head())

print(adata_sub.var_names[:10])

adata_sub.write_h5ad("../../data/for_saigeQTL/untreated_subset_for_saigeQTL.h5ad")
### getting the 
# geno_pheno_mapping.tsv

# get unique donor IDs
donors = adata_sub.obs["donor_id"].unique()

df = pd.DataFrame({
    "Genotype": donors,
    "RNA": donors
})

print(df.head())
print(df["Genotype"].nunique())
#242
print(df.shape)

df.to_csv(
    "/lustre/scratch125/humgen/projects_v2/otar2065/OTAR2065_sc_eQTL/data/for_saigeQTL/geno_pheno_mapping.tsv",
    sep="\t",
    index=False
)
