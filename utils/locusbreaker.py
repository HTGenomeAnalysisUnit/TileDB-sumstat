import pandas as pd
from typing import List


def locus_breaker(
    tiledb_data,
    pvalue_sig: float = 5e-8,
    pvalue_limit: float = 5e-6,
    hole_size: int = 250000
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
    'cell_type': pd.Series(dtype='object'),
    'position': pd.Series(dtype='int32'),
    'gene': pd.Series(dtype='object'),
    'SNP': pd.Series(dtype='object'),
    'af': pd.Series(dtype='float32'),
    'se': pd.Series(dtype='float32'),
    'beta': pd.Series(dtype='float32'),
    'p-value': pd.Series(dtype='float64'),
    'allele0': pd.Series(dtype='object'),
    'allele1': pd.Series(dtype='object')
}

    # Filter rows based on the p_limit threshold
    tiledb_data = tiledb_data[tiledb_data["p-value"] < pvalue_limit]
    
    # If no rows remain after filtering, return an empty DataFrame
    if tiledb_data.empty:
        return pd.DataFrame(expected_schema)
    # Group by gene first, then calculate regions within each gene
    trait_res = []
    trait_res_allsnp = []
    for gene, gene_df in tiledb_data.groupby("gene"):
        # Find regions where gaps between positions exceed hole_size within each chromosome
        gaps = gene_df["position"].diff() > hole_size
        group = gaps.cumsum()
        # Group by the identified regions within the chromosome
        for _, group_df in gene_df.groupby(group):
            if group_df["p-value"].min() < pvalue_sig:
                start_pos = group_df["position"].min()
                end_pos = group_df["position"].max()
                best_snp = group_df.loc[group_df["p-value"].idxmin()]

                # Store the interval with the best SNP
                line_res = [gene, start_pos, end_pos, best_snp["position"], best_snp["p-value"]] + best_snp.tolist()
                trait_res.append(line_res)

                #Collect all SNPs within the region
                for _, snp_row in group_df.iterrows():
                    snp_res = [gene, snp_row["position"], snp_row["p-value"]] + snp_row.tolist()
                    trait_res_allsnp.append(snp_res)

    # Convert results to a DataFrame

    columns = ["gene", "start", "end", "snp_pos", "snp_pval"] + tiledb_data.columns.tolist()
    trait_res_df = pd.DataFrame(trait_res, columns=columns)
    trait_res_df = trait_res_df.drop(trait_res_df.columns[0], axis=1)
    trait_res_df = trait_res_df.drop(columns=["position", "p-value"])
        
    columns = ["gene", "snp_pos", "snp_pval"] + tiledb_data.columns.tolist()
    trait_res_allsnp_df = pd.DataFrame(trait_res_allsnp, columns=columns)
    #trait_res_allsnp_df = trait_res_allsnp_df.drop(trait_res_allsnp_df.columns[0], axis=1)
    trait_res_allsnp_df = trait_res_allsnp_df.drop(columns=["snp_pos", "snp_pval"])

    return [trait_res_df,trait_res_allsnp_df]