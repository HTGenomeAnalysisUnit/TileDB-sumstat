"""Locusbreaker-based export from TileDB."""
import os
import random

import pandas as pd
import polars as pl
import tiledb

from tdbsumstat.utils.locusbreaker_plpl import locusbreaker_plpl


def _query_tiledb_for_locusbreaker(
    uri_path: str,
    chrom: int,
    type_sumstat: str,
    trait: str = None,
    cell: str = None,
    gene: str = None,
) -> "pa.Table":
    """Query TileDB for a specific chromosome/trait combination.

    Returns a PyArrow table for use with locusbreaker_plpl.
    """
    with tiledb.open(uri_path, mode="r") as tiledb_data:
        if type_sumstat == "gwas":
            return tiledb_data.query(dims=["CHR", "TRAIT", "POS"]).df[chrom, trait, :]
        else:
            return tiledb_data.query(dims=["CHR", "CELL", "GENE", "POS"], return_arrow=True).df[
                chrom, cell, gene, :
            ]


def export_with_locusbreaker(
    uri_path: str,
    df_meta: pl.DataFrame,
    table_lb: str,
    maf_lb: float,
    hole_lb: int,
    locus_max_size_lb: int,
    cis_trans_lb: str,
    type_sumstat: str,
    out: str,
    batch_name: str,
    pvalue_sig: float = 5e-8,
    pvalue_limit: float = 5e-6,
) -> None:
    """Run locusbreaker on TileDB data and export loci intervals and segments.

    Parameters
    ----------
    uri_path : str
        Path to the TileDB array.
    df_meta : pl.DataFrame
        Polars DataFrame with merged metadata (used by locusbreaker_plpl).
    table_lb : str
        Path to a CSV table with CHR, TRAIT (and optionally SIG, LIM columns).
    maf_lb : float
        MAF filter applied before locusbreaker.
    hole_lb : int
        Minimum base-pair distance to separate loci.
    locus_max_size_lb : int
        Maximum allowed locus size in base pairs.
    cis_trans_lb : str
        Filter type: "cis" or "trans" (for QTL data).
    type_sumstat : str
        Type of summary statistics: "gwas" or "qtl".
    out : str
        Output file prefix.
    batch_name : str
        Batch identifier appended to output file names.
    pvalue_sig : float
        P-value threshold for significant SNPs (default: 5e-8).
    pvalue_limit : float
        P-value threshold for locus boundary definition (default: 5e-6).
    """
    print("Starting LocusBreaker")
    traits = pd.read_csv(table_lb)
    traits = traits.astype({"CHR": "int16"})

    if not batch_name:
        batch_name = random.randint(1, 10000000)

    for _index, trait in traits.iterrows():
        if "SIG" in traits.columns:
            pvalue_sig = trait["SIG"]
            pvalue_limit = trait["LIM"]

        if type_sumstat == "gwas":
            query = _query_tiledb_for_locusbreaker(uri_path, trait["CHR"], type_sumstat, trait=trait["TRAIT"])
        else:
            cell, genes = trait["TRAIT"].split(";")
            query = _query_tiledb_for_locusbreaker(uri_path, trait["CHR"], type_sumstat, cell=cell, gene=genes)

        result = locusbreaker_plpl(
            query,
            maf=maf_lb,
            pvalue_sig=pvalue_sig,
            pvalue_limit=pvalue_limit,
            locus_max_size=locus_max_size_lb,
            hole_size=hole_lb,
            cis_trans_lb=cis_trans_lb,
            type_sumstat=type_sumstat,
            metadata=df_meta,
        )

        if not len(result) == 0 and not result[0].empty:
            if result and isinstance(result[0], pd.DataFrame) and not result[0].shape[0] == 0:
                interval = result[0]
                segments = result[1]

                write_header_interval = not os.path.exists(f"{out}_batch_{batch_name}_interval.csv")
                write_header_segment = not os.path.exists(f"{out}_batch_{batch_name}_segment.csv")
                interval.to_csv(
                    f"{out}_batch_{batch_name}_interval.csv",
                    mode="a",
                    index=False,
                    header=write_header_interval,
                )
                segments.to_csv(
                    f"{out}_batch_{batch_name}_segment.csv",
                    mode="a",
                    index=False,
                    header=write_header_segment,
                )
