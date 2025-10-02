import pandas as pd
import tiledb
rows = []
tdb = tiledb.open('/lustre/scratch124/humgen/projects_v2/cardinal_analysis/analysis/core_dataset/freeze3/tiledbs/Tensor/TileDB_f3_UKB_cis_t2', 'w')
# Loop over all celltypes in metadata
for celltype in merged_metadata["celltype"]:
    if celltype in merged_metadata:  # make sure the key exists in dict
        groups = merged_metadata[celltype]
        for group, values in groups.items():  # e.g. group "11"
            for gene_id in values:
                rows.append({
                    "CHR": group,
                    "CELL": celltype,
                    "GENE": gene_id,
                    "N": values[gene_id]["N"],
                    "ACAT": values[gene_id]["ACAT"],
                    "PHENOVAR": values[gene_id]["PHENO_VAR"],
                })
df = pd.DataFrame(rows)

tdb.meta["metadata"] = json.dumps(merged_metadata)
df.to_csv("/lustre/scratch124/humgen/projects_v2/cardinal_analysis/analysis/core_dataset/freeze3/tiledbs/Tensor/Bangladeshi/metadata_gh_bangladeshi_f3_celltype1_acat.csv",index = False)