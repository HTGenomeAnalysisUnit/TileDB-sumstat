import json
import re
from collections import defaultdict

# === STEP 1: Read the entire string from file ===
with open("/project/cardinal/QTLs/freeze3/TileDBs/TileDB_metanalyses_gh_celltype2_freeze3_metadata.json", "r", encoding="utf-8") as f:
    data_str = f.read()

# === STEP 2: Split based on the start of each JSON fragment ===
# We capture the first '{' before each file_path
fragments = re.findall(r'(\{[^{]*"file_path":.*?)(?=\{[^{]*"file_path":|\Z)', data_str, flags=re.S)

print(f"🧩 Found {len(fragments)} fragments")

# === STEP 3: Parse each fragment ===
parsed_fragments = []
for i, frag in enumerate(fragments):
    try:
        parsed_fragments.append(json.loads(frag))
    except json.JSONDecodeError as e:
        print(f"⚠️ Fragment {i} failed to parse at char {e.pos}: {e}")
        print(frag[:300] + "...")
        # optional: skip or raise
        continue

# === STEP 4: Merge fragments as before ===
merged = {
    "trait": [],
    "celltype": [],
}
celltype_data = defaultdict(lambda: defaultdict(dict))

for frag in parsed_fragments:
    if "trait" in frag and isinstance(frag["trait"], list):
        merged["trait"].extend(frag["trait"])
    if "celltype" in frag and isinstance(frag["celltype"], list):
        merged["celltype"].extend(frag["celltype"])
    for key, val in frag.items():
        if key in ("trait", "celltype", "CELL", "file_path"):
            continue
        for num, genes in val.items():
            for gene_id, gene_data in genes.items():
                celltype_data[key][num][gene_id] = gene_data

for celltype, numbers in celltype_data.items():
    merged[celltype] = numbers

merged["trait"] = list(dict.fromkeys(merged["trait"]))
merged["celltype"] = list(dict.fromkeys(merged["celltype"]))

# === STEP 5: Save merged JSON ===
with open("metadata_merged.json", "w", encoding="utf-8") as out_f:
    json.dump(merged, out_f, indent=2)

print("✅ Merged JSON written to metadata_merged.json")
