"""Metadata export and recomputation from TileDB."""
import json

import pandas as pd
import polars as pl
import tiledb

from tdbsumstat.utils import acat_optimized


def export_metadata(uri_path: str, type_sumstat: str, out: str) -> None:
    """Export merged metadata from TileDB to a CSV file.

    Parameters
    ----------
    uri_path : str
        Path to the TileDB array.
    type_sumstat : str
        Type of summary statistics: "gwas" or "qtl".
    out : str
        Output file prefix (``{out}_meta.csv`` will be created).
    """
    tiledb_db = tiledb.open(uri_path, mode="r")
    tiledb_meta = json.loads(tiledb_db.meta["merged_metadata"])
    rows = []

    if type_sumstat == "qtl":
        for cell_type in tiledb_meta["CELL"]:
            samples = tiledb_meta.get(cell_type, {})
            for sample_id, genes in samples.items():
                for gene_id, metrics in genes.items():
                    rows.append({
                        "CELL": cell_type,
                        "CHR": sample_id,
                        "GENE": gene_id,
                        **metrics,
                    })
    else:
        for trait in tiledb_meta["traits"]:
            samples = tiledb_meta.get(trait, {})
            for chrom, genes in samples.items():
                for gene_id, metrics in genes.items():
                    rows.append({
                        "TRAIT": trait,
                        "CHR": chrom,
                        **metrics,
                    })

    df = pd.DataFrame(rows)
    df.to_csv(f"{out}_meta.csv", index=False)
    tiledb_db.close()


def recompute_metadata(
    tiledb_export: tiledb.Array,
    df_meta: pl.DataFrame,
    trait_list: str,
    type_sumstat: str,
    mac: int,
    out: str,
    batch_name: str,
) -> None:
    """Recompute metadata statistics after applying MAC filter and save to CSV.

    This does **not** modify data stored inside the TileDB array.

    Parameters
    ----------
    tiledb_export : tiledb.Array
        Open TileDB array in read mode.
    df_meta : pl.DataFrame
        Polars DataFrame with merged metadata.
    trait_list : str
        Path to a CSV file with traits/cells to recompute.
    type_sumstat : str
        Type of summary statistics: "gwas" or "qtl".
    mac : int
        Minimum minor allele count filter.
    out : str
        Output file prefix.
    batch_name : str
        Batch identifier appended to output file names.
    """
    trait_list_pd = pd.read_csv(trait_list)

    for _record, trait in trait_list_pd.iterrows():
        if type_sumstat == "gwas":
            tiledb_query = tiledb_export.query().df[int(trait["CHR"]), trait["TRAIT"].to_string(), :]
            n = df_meta.filter(pl.col("TRAIT") == trait["TRAIT"]).select("N")["N"][0]
        else:
            print(trait["CELL"])
            tiledb_query = tiledb_export.query().df[int(trait["CHR"]), trait["CELL"], :, :]
            n = df_meta.filter(pl.col("CELL") == trait["CELL"]).select("N")["N"][0]
            print(n)

        if "N" not in tiledb_query.columns:
            tiledb_query["N"] = n
            print(tiledb_query)

        tiledb_query_pl = pl.from_pandas(tiledb_query)
        tiledb_query_pl = tiledb_query_pl.with_columns(
            (2 * pl.col("N") * pl.min_horizontal("EAF", (1 - pl.col("EAF")))).alias("MAC")
        ).filter(pl.col("MAC") > mac)

        if type_sumstat == "gwas":
            chr_gene_agg = tiledb_query_pl.group_by(["CHR", "TRAIT"]).agg([
                pl.col("P")
                .map_batches(
                    lambda s: pl.Series([acat_optimized(s)]),
                    return_dtype=pl.Float64,
                )
                .alias("ACAT_LIST"),
                pl.col("N").first().alias("N"),
                pl.min("P").alias("MIN_P"),
            ])
            chr_gene_agg = chr_gene_agg.with_columns(
                pl.col("ACAT_LIST").list.first().alias("ACAT")
            )
        else:
            chr_gene_agg = tiledb_query_pl.group_by(["CHR", "CELL", "GENE"]).agg([
                pl.col("P")
                .map_batches(
                    lambda s: pl.Series([acat_optimized(s)]),
                    return_dtype=pl.Float64,
                )
                .alias("ACAT_LIST"),
                pl.col("N").first().alias("N"),
                pl.min("P").alias("MIN_P"),
            ])
            chr_gene_agg = chr_gene_agg.with_columns(
                pl.col("ACAT_LIST").list.first().alias("ACAT"),
            ).drop("ACAT_LIST")

        chr_gene_agg_pd = chr_gene_agg.to_pandas()
        chr_gene_agg_pd.to_csv(f"{out}_batch_{batch_name}_metadata.csv", mode="a", index=False)
