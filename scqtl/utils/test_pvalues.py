import tiledb 
import pandas as pd
import polars as pl
tdb = tiledb.open("/project/cardinal/QTLs/freeze3/TileDBs/mashr/TileDB_mashr_ukb_celltype2_freeze3_large_chunks")
metadata = pd.read_csv("/group/soranzo/bruno.ariano/ingest_mashr_f3/mashr_metadata_missing.csv")
metadata_prior = pd.read_csv("/group/soranzo/bruno.ariano/ingest_mashr_f3/mashr_metadata_pvalues.csv")

for record_index, record in df.iterrows():
    chr = record["CHR"]
    cell,gene = record["CELL:GENE"].split(':')
    df_query = tdb.query(dims=["GENE"], attrs = ["P"]).df_query[chr,cell,gene]
    min_p = min(df_query["P"])
    print(chr)
    print(cell)
    print(min_p)

import numpy as np

# Create a new column in metadata to store min P values
metadata_df["MIN_P"] = np.nan

for record_index, record in metadata_df.iterrows():
        chr = int(record["CHR"])
        cell, gene = record["CELL:GENE"].split(':')
        # Query the TileDB array for the given coordinates
        df_query = tdb.query(dims=["GENE"], attrs = ["P"]).df[chr,cell,gene,:]
        # Take the minimum p-value
        min_p = df_query["P"].min()
        metadata_df.at[record_index, "MIN_P"] = min_p
        print(chr)
        print(cell)
        print(min_p)
pd.concat[metadata_prior,metadata]




metadata_df["MIN_P"] = np.nan
metadata_df[["CELL", "GENE"]] = metadata_df["CELL:GENE"].str.split(":", expand=True)

groups = list(metadata_df.groupby(["CHR", "CELL"]))  # materialize groups
total_groups = len(groups)

for i, ((chr_, cell), group) in enumerate(groups, start=1):
    print(f"Processing group {i}/{total_groups} → CHR={chr_}, CELL={cell}", flush=True)
    genes = group["GENE"].tolist()
    df_query = tdb.query(dims=["GENE"], attrs=["P"]).df[int(chr_), cell, genes, :]
    min_p_map = df_query.groupby("GENE")["P"].min().to_dict()
    metadata_df.loc[group.index, "MIN_P"] = group["GENE"].map(min_p_map)


all_min_chunks = []

for chr_ in metadata_df["CHR"].unique():
        print(f"Reading CHR {chr_}...")
        df_chr = tdb.query(attrs=["P"], dims=["CHR", "CELL", "GENE"]).df[int(chr_), :, :, :]
        df_min_chr = df_chr.groupby(["CHR", "CELL", "GENE"], as_index=False)["P"].min()
        all_min_chunks.append(df_min_chr)

metadata_df[["CELL", "GENE"]] = metadata_df["CELL:GENE"].str.split(":", expand=True)
metadata_df["CHR"] = metadata_df["CHR"].astype(int)

metadata_df = metadata_df.merge(
    df_min,
    on=["CHR", "CELL", "GENE"],
    how="left",          # keep all rows from merged_metadata
    suffixes=("", "_min")
)




metadata_df_complete[["CELL", "GENE"]] = metadata_df_complete["CELL:GENE"].str.split(":", expand=True)

for col in ["CHR", "CELL", "GENE"]:
    metadata_df[col] = metadata_df[col].astype(str)
    metadata_df_complete[col] = metadata_df_complete[col].astype(str)

# Merge only the MIN_P column from df_with_minp
merged = metadata_df_complete.merge(
    metadata_df[["CHR", "CELL", "GENE", "MIN_P"]],
    on=["CHR", "CELL", "GENE"],
    how="left",
    suffixes=("", "_new")
)

# Fill missing MIN_P values with the ones from df_with_minp
merged["MIN_P"] = merged["MIN_P"].fillna(merged["MIN_P_new"])

# Drop the temporary column
merged.drop(columns=["MIN_P_new"], inplace=True)



with open("/project/cardinal/QTLs/freeze3/TileDBs/mashr/TileDB_mashr_ukb_celltype2_freeze3_large_chunks_metadata.json", "r") as f:
    data = json.load(f)





import json

def parse_multiple_json(s):
    objects = []
    decoder = json.JSONDecoder()
    while s:
        s = s.lstrip()  # Skip leading whitespace
        if not s:
            break
        try:
            obj, idx = decoder.raw_decode(s)
            objects.append(obj)
            s = s[idx:]
        except json.JSONDecodeError as e:
            print(f"Stopped parsing at position {e.pos}: {e}")
            print(f"Remaining unparsed data starts with: {s[:200]}...")
            break
    return objects

# Usage in your REPL
file_path = "/lustre/scratch124/humgen/projects_v2/cardinal_analysis/analysis/core_dataset/freeze3/tiledbs/Tensor/TileDB_mashr_ukb_celltype2_freeze3_large_chunks_metadata.json"

with open(file_path, "r") as f:
    content = f.read()

all_objects = parse_multiple_json(content)
print(f"Parsed {len(all_objects)} JSON objects.")
# all_objects[0] is the main metadata, all_objects[1:] are the gene-specific data