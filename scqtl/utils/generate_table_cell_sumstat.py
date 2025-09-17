import pandas as pd
import tiledb
rows = []
tdb = tiledb.open('/lustre/scratch124/humgen/projects_v2/cardinal_analysis/analysis/core_dataset/freeze3/tiledbs/Tensor/TileDB_f3_Bangladeshi_cis_celltype1', 'w')
# Loop over all celltypes in metadata
for celltype in merged_metadata["celltype"]:
    if celltype in merged_metadata:  # make sure the key exists in dict
        groups = merged_metadata[celltype]
        for group, values in groups.items():  # e.g. group "11"
            for gene_id, value in values:
                rows.append({
                    "CHR": group,
                    "CELL": celltype,
                    "GENE": gene_id,
                    "ACAT": value
                })
df = pd.DataFrame(rows)

tdb.meta["metadata"] = json.dumps(merged_metadata)
df.to_csv("/lustre/scratch124/humgen/projects_v2/cardinal_analysis/analysis/core_dataset/freeze3/tiledbs/Tensor/metadata_gh_bangladeshi_f3_celltype1_acat.csv",index = False)