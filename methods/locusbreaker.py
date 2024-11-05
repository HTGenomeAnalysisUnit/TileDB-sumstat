from dask import delayed, compute
import pandas as pd
from typing import List


@delayed
def locus_breaker(
    tiledb_data,
    out:  str,
    pvalue_sig: float = 5e-8,
    pvalue_limit: float = 5e-6,
    hole_size: int = 250000,
    column_list_select: List[str] = [
        "position",
        "SNP_ID",
        "BETA",
        "p-value",
    ],
    map_attributes: dict = None
) -> pd.DataFrame:
    """
    Breaking genome in locus
    Returns a series of parquet files describing the loci created dynamically around significant SNPs.
    :param tiledb_data: TileDBVCF data (default: None)
    :param pvalue_sig: P-value threshold in -log10 format used to create the regions around significant SNPs (default: 5)
    :param pvalue_limit: P-value threshold in -log10 format for loci borders (default: 5)
    :param hole_size: Minimum pair-base distance between SNPs in different loci (default: 250000)
    :return: DataFrame with the loci information
    """
    expected_schema = {
    'snp_pos': pd.Series(dtype='int64'),
    'snp_PVAL': pd.Series(dtype='float64'),
    'beta': pd.Series(dtype='object'),
    'p-value': pd.Series(dtype='object'),
    'SNP': pd.Series(dtype='int64')
}

    # Convert fmt_LP from list to float
    if tiledb_data.empty:
        print("this region is empty")
        return pd.DataFrame(expected_schema)
    
    #tiledb_data["PVAL"] = tiledb_data["PVAL"].apply(lambda x: float(x[0]))

    # Filter rows based on the p_limit threshold
    tiledb_data = tiledb_data[tiledb_data["p-value"] < pvalue_limit]

    # If no rows remain after filtering, return an empty DataFrame
    if tiledb_data.empty:
        return pd.DataFrame(expected_schema)

    # Group by 'contig' (chromosome) first, then calculate regions within each chromosome
    trait_res = []

    gaps = tiledb_data["position"].diff() > hole_size
    group = gaps.cumsum()

    # Group by the identified regions within the chromosome
    for _, group_df in tiledb_data.groupby(group):
        if group_df["p-value"].min() < pvalue_sig:
            start_pos = group_df["position"].min()
            end_pos = group_df["position"].max()
            best_snp = group_df.loc[group_df["p-value"].idxmin()]

            # Store the interval with the best SNP
            #line_res = [start_pos, end_pos, best_snp["position"], best_snp["p-value"]] + best_snp.tolist()
            #trait_res.append(line_res)

            # Collect all SNPs within the region
            for _, snp_row in group_df.iterrows():
                snp_res = [start_pos, end_pos, snp_row["position"], snp_row["p-value"]] + snp_row.tolist()
                trait_res.append(snp_res)

    # Convert results to a DataFrame
    columns = ["start", "end", "snp_pos", "snp_PVAL"] + tiledb_data.columns.tolist()
    trait_res_df = pd.DataFrame(trait_res, columns=columns)

    # Drop specific columns including 'start' and 'end'
    #trait_res_df = trait_res_df.drop(columns=["position", "p-value", "start", "end"])

    # Remove one of the duplicate 'contig' columns if present
    trait_res_df = trait_res_df.loc[:, ~trait_res_df.columns.duplicated()]
    trait_res_df.to_parquet(out, engine = "pyarrow")

    #columns_attribute_mapping = {v: k for k, v in map_attributes.items() if v in trait_res_df.columns}

    #trait_res_df.rename(columns=columns_attribute_mapping)

    # Rename the columns using the map

    return trait_res_df
