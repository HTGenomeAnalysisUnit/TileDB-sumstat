"""SNP-based export from TileDB."""
import numpy as np
import pandas as pd
import polars as pl
import tiledb


def export_by_snp(
    tiledb_export: tiledb.Array,
    snp: str,
    attr: str,
    type_sumstat: str,
    out: str,
) -> None:
    """Query TileDB by a list of SNPs and export results to CSV.

    Parameters
    ----------
    tiledb_export : tiledb.Array
        Open TileDB array in read mode.
    snp : str
        Path to a CSV file with columns CHR, POS, TRAIT (and optionally GENE for QTL).
    attr : str
        Comma-separated list of attributes to export.
    type_sumstat : str
        Type of summary statistics: "gwas" or "qtl".
    out : str
        Output file prefix.
    """
    snp_list = pd.read_csv(snp, dtype={"CHR": int, "POS": np.uint32, "TRAIT": str})

    if type_sumstat == "gwas":
        trait_list = snp_list["TRAIT"].unique().tolist()
    else:
        snp_list[["CELL", "GENE"]] = snp_list["TRAIT"].str.split(":", n=2, expand=True)
        trait_list = snp_list["CELL"].unique().tolist()

    for trait in trait_list:
        if type_sumstat == "gwas":
            chrom_list = snp_list[snp_list["TRAIT"] == trait]["CHR"].unique().tolist()
        else:
            chrom_list = snp_list[snp_list["CELL"] == trait]["CHR"].unique().tolist()

        for chrom in chrom_list:
            if type_sumstat == "gwas":
                snp_list_refined = (
                    snp_list[(snp_list["CHR"] == chrom) & (snp_list["TRAIT"] == trait)]["POS"]
                    .unique()
                    .tolist()
                )
                tiledb_query = tiledb_export.query(
                    attrs=attr.split(","),
                    return_arrow=True,
                ).df[chrom, trait, :]
                tiledb_query_pl = pl.from_arrow(tiledb_query)
                tiledb_query_pd = tiledb_query_pl.filter(pl.col("POS").is_in(snp_list_refined)).to_pandas()
            else:
                snp_list_refined = (
                    snp_list[(snp_list["CHR"] == chrom) & (snp_list["CELL"] == trait)]["POS"]
                    .unique()
                    .tolist()
                )
                gene_list = (
                    snp_list[(snp_list["CHR"] == chrom) & (snp_list["CELL"] == trait)]["GENE"]
                    .unique()
                    .tolist()
                )
                tiledb_query = tiledb_export.query(attrs=attr.split(",")).df[chrom, trait, gene_list, :]
                tiledb_query_pl = pl.from_arrow(tiledb_query)
                tiledb_query_pd = tiledb_query_pl.filter(pl.col("POS").is_in(snp_list_refined)).to_pandas()

            tiledb_query_pd.to_csv(f"{out}_{trait}_{chrom}.csv", mode="a", index=False)
