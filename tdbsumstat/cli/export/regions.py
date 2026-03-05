"""Region-based export from TileDB."""
import pandas as pd
import tiledb


def export_by_regions(
    tiledb_export: tiledb.Array,
    table_regions: str,
    attr: str,
    type_sumstat: str,
    out: str,
) -> None:
    """Query TileDB by genomic regions table and export results to CSV.

    Parameters
    ----------
    tiledb_export : tiledb.Array
        Open TileDB array in read mode.
    table_regions : str
        Path to a CSV file with columns CHR, START, END, TRAIT.
    attr : str
        Comma-separated list of attributes to export.
    type_sumstat : str
        Type of summary statistics: "gwas" or "qtl".
    out : str
        Output file path.
    """
    pd_region = pd.read_csv(table_regions)
    counter_nonempty_region = 0

    for _ind, row in pd_region.iterrows():
        if type_sumstat == "gwas":
            trait = row["TRAIT"]
            region = tiledb_export.query(
                dims=["CHR", "POS", "TRAIT"],
                attrs=attr.split(","),
            ).df[int(row["CHR"]), trait, int(row["START"]):int(row["END"])]
        else:
            cell, gene = row["TRAIT"].split(":")
            region = tiledb_export.query(
                dims=["CHR", "POS", "CELL", "GENE"],
                attrs=attr.split(","),
            ).df[int(row["CHR"]), cell, gene, int(row["START"]):int(row["END"])]

        if len(region) > 0:
            if counter_nonempty_region == 0:
                region.to_csv(out, mode="a", index=False, header=True)
                counter_nonempty_region += 1
            else:
                region.to_csv(out, mode="a", index=False, header=False)
