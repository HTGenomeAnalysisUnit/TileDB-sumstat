"""Common TileDB helper functions for export operations."""
import json

import polars as pl
import tiledb


def open_tiledb_and_load_metadata(uri_path: str, type_sumstat: str):
    """Open TileDB array and load metadata as a Polars DataFrame.

    Parameters
    ----------
    uri_path : str
        Path to the TileDB array.
    type_sumstat : str
        Type of summary statistics: "gwas" or "qtl".

    Returns
    -------
    tuple[tiledb.Array, pl.DataFrame]
        Open TileDB array (read mode) and a Polars DataFrame with metadata.
    """
    tiledb_export = tiledb.open(uri_path, mode="r")
    metadata = json.loads(tiledb_export.meta["merged_metadata"])
    rows = []
    if type_sumstat == "qtl":
        for cell in metadata["CELL"]:
            for chrom, genes in metadata[cell].items():
                for gene, stats in genes.items():
                    rows.append({
                        "CELL": cell,
                        "CHR": chrom,
                        "GENE": gene,
                        **stats
                    })
        df_meta = pl.DataFrame(rows)
        df_meta = df_meta.with_columns(pl.col("CHR").cast(pl.UInt16))
    else:
        for trait in metadata["trait"]:
            rows.append({
                "TRAIT": trait,
                **metadata[trait]
            })
        df_meta = pl.DataFrame(rows)
    return tiledb_export, df_meta
