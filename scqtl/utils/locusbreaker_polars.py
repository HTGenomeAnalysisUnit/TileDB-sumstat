import polars as pl
from typing import List
from utils import compute_pheno_variance
import pandas as pd
def locus_breaker(
    tiledb_data: pl.DataFrame,
    pvalue_sig: float = 5e-8,
    pvalue_limit: float = 5e-6,
    hole_size: int = 250000,
    phenovar: bool = False
) -> List[pl.DataFrame]:
    """
    Breaking genome into loci.
    Returns a list of DataFrames describing the loci created dynamically around significant SNPs.

    :param tiledb_data: TileDBVCF data as a Polars DataFrame.
    :param pvalue_sig: P-value threshold to define significant SNPs.
    :param pvalue_limit: P-value threshold for loci borders.
    :param hole_size: Minimum base-pair distance between SNPs in different loci.
    :param phenovar: Boolean indicating whether to compute phenotype variance.
    :return: List of DataFrames with the loci information.
    """

    # Convert to Polars DataFrame if not already
    if not isinstance(tiledb_data, pl.DataFrame):
        tiledb_data_pl = pl.from_pandas(tiledb_data)

    # Filter by p-value limit
    tiledb_data_pl = tiledb_data_pl.filter(pl.col("P") < pvalue_limit)

    # Compute phenotype variance if requested
    tiledb_data_pl = tiledb_data_pl.with_columns(
        pl.lit(1.0).alias("S") if not phenovar else pl.Series("S", compute_pheno_variance(tiledb_data_pl))
    )

    # If no rows remain, return empty DataFrames
    if tiledb_data_pl.is_empty():
        return [pl.DataFrame([]), pl.DataFrame([])]

    trait_res = []
    trait_res_allsnp = []

    for gene, gene_df in tiledb_data_pl.group_by("GENE"):
        gene_df = gene_df.sort("POS")  # Ensure sorted by position

        # Identify loci based on gap size
        gaps = gene_df["POS"].diff().fill_null(hole_size + 1) > hole_size
        group = gaps.cum_sum()

        # Group SNPs into loci
        for _, group_df in gene_df.with_columns(group.alias("group")).group_by("group"):
            if group_df["P"].min() < pvalue_sig:
                start_pos = group_df["POS"].min()
                end_pos = group_df["POS"].max()
                best_snp = group_df.sort("P").head(1)

                # Extract values while ensuring correct schema
                best_snp_values = best_snp.row(0)[1:]  # Drop the first tuple element
                line_res = [gene[0], start_pos, end_pos, best_snp["POS"][0], best_snp["P"][0]] + list(best_snp_values)

                trait_res.append(line_res)

                # Collect all SNPs from this locus
                for row in group_df.iter_rows():
                    snp_res = [gene, row[1], row[2]] + list(row[0:-1])
                    trait_res_allsnp.append(snp_res)

    # Define column names
    columns = ["GENE", "start", "end", "snp_pos", "snp_pval"] + tiledb_data.columns.tolist() + ["S"]
    trait_res_df = pd.DataFrame(trait_res, columns=columns)
    columns = ["GENE", "snp_pos", "snp_pval"] + tiledb_data.columns.tolist() + ["S"]
    trait_res_allsnp_df = pd.DataFrame(trait_res_allsnp, columns=columns)
    trait_res_allsnp_df = trait_res_allsnp_df.drop(columns=["snp_pos", "snp_pval"])

    # Create DataFrames while ensuring the correct structure
    # Ensure "GENE" is not duplicated in the schema
    trait_res_df = trait_res_df.drop(trait_res_df.columns[0], axis=1)
    trait_res_df = trait_res_df.drop(columns=["POS", "P"])

    #trait_res_allsnp_df = trait_res_allsnp_df.drop(trait_res_allsnp_df.columns[0], axis=1)
    return [trait_res_df, trait_res_allsnp_df]
