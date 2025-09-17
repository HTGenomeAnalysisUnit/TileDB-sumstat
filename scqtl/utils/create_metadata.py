import json
from collections import defaultdict
input_path = "/lustre/scratch124/humgen/projects_v2/cardinal_analysis/analysis/core_dataset/freeze3/tiledbs/Tensor/TileDB_f3_Bangladeshi_cis_celltype1_metadata.json"
output_path = input_path.replace(".json", "_fixed.json")

objs = []
with open(input_path) as f:
    text = f.read()

decoder = json.JSONDecoder()
idx = 0
while idx < len(text):
    obj, idx = decoder.raw_decode(text, idx)
    objs.append(obj)

with open(output_path, "w") as f:
    json.dump(objs, f, indent=2)

print(f"✅ Salvaged {len(objs)} metadata records into {output_path}")


fixed_path = "/lustre/scratch124/humgen/projects_v2/cardinal_analysis/analysis/core_dataset/freeze3/tiledbs/Tensor/TileDB_f3_Bangladeshi_cis_celltype1_metadata_fixed.json"
output_path = fixed_path.replace("_fixed.json", "_merged.json")

with open(fixed_path) as f:
    metadata_records = json.load(f)

# Create a merged dictionary
merged_metadata = {
    "file_path": [],
    "trait": [],
    "celltype": []
}

# For storing per-cell/chromosome info if present
cell_chrom_data = defaultdict(lambda: defaultdict(list))

for record in metadata_records:
    # Merge file paths
    if "file_path" in record and record["file_path"] not in merged_metadata["file_path"]:
        merged_metadata["file_path"].append(record["file_path"])
    
    # Merge traits
    if "trait" in record:
        for t in record["trait"]:
            if t not in merged_metadata["trait"]:
                merged_metadata["trait"].append(t)
    
    # Merge celltypes and their per-chromosome data
    if "celltype" in record:
        for cell in record["celltype"]:
            if cell not in merged_metadata["celltype"]:
                merged_metadata["celltype"].append(cell)
            # Merge chromosome-level info if exists
            if cell in record:
                for chrom, genes in record[cell].items():
                    existing_genes = {tuple(g) for g in cell_chrom_data[cell][chrom]}
                    for g in genes:
                        if tuple(g) not in existing_genes:
                            cell_chrom_data[cell][chrom].append(g)

# Add per-cell/chromosome info
for cell, chrom_dict in cell_chrom_data.items():
    merged_metadata[cell] = chrom_dict

# Save as a single JSON object
with open(output_path, "w") as f:
    json.dump(merged_metadata, f, indent=2)

print(f"✅ Merged metadata saved to {output_path}")
