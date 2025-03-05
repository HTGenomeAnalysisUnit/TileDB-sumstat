import pandas as pd
from typing import List
from utils import compute_pheno_variance

def locus_breaker(
    tiledb_data,
    pvalue_sig: float = 5e-8,
    pvalue_limit: float = 5e-6,
    hole_size: int = 250000,
    phenovar: bool = False,
    maf: float = 0.01,
    category: bool = False
) -> pd.DataFrame:
    """
    Breaking genome in loci and returning all SNPs within the defined loci boundaries.
    :param tiledb_data: TileDBVCF data
    :param pvalue_sig: P-value threshold for significant SNPs
    :param pvalue_limit: P-value threshold for defining loci
    :param hole_size: Minimum base-pair distance to separate loci
    :param phenovar: Compute phenotypic variance or not
    :param maf: Minor allele frequency threshold
    :param category: cis/trans filter for SNPs based on distance
    :param expansion_size: Number of base pairs to expand loci boundaries
    :return: Two DataFrames, one with loci regions and another with all SNPs in loci
    """
    # Create a copy of the original dataset before filtering
    original_data = tiledb_data.copy()
    
    if phenovar:
        original_data["S"] = compute_pheno_variance(original_data)
    else:
        original_data["S"] = 1.0
    
    # Filter for SNPs below the p-value limit to define loci
    loci_snps = tiledb_data[tiledb_data["P"] < pvalue_limit].copy()
    
    if loci_snps.empty:
        return []
    
    # Apply cis/trans filtering if needed
    if category == "cis":
        loci_snps = loci_snps[(loci_snps["DIST"] > -1000000) & (loci_snps["DIST"] < 1000000)]
    elif category == "trans":
        loci_snps = loci_snps[(loci_snps["DIST"] < -1000000) | (loci_snps["DIST"] > 1000000)]
    
    trait_res = []
    all_snp_res = []
    
    for gene, gene_df in loci_snps.groupby("GENE"):
        gaps = gene_df["POS"].diff() > hole_size
        group = gaps.cumsum()
        
        for _, group_df in gene_df.groupby(group):
            if group_df["P"].min() < pvalue_sig:
                start_pos = group_df["POS"].min() - 100000
                end_pos = group_df["POS"].max() + 100000
                best_snp = group_df.loc[group_df["P"].idxmin()]
                region = f"{group_df['CHR'].iloc[0]}:{start_pos}:{end_pos}"
                trait_id = f"{group_df['CELL']}_{gene}"
                
                trait_res.append([trait_id, start_pos, end_pos, best_snp["POS"], best_snp["P"]] + best_snp.tolist())
                
                # Include all SNPs within the expanded region from the original dataset
                expanded_snps = original_data[
                    (original_data["CHR"] == group_df["CHR"].iloc[0]) &
                    (original_data["POS"] >= start_pos) &
                    (original_data["POS"] <= end_pos)
                ]
                for _, snp_row in expanded_snps.iterrows():
                    all_snp_res.append([trait_id, region, snp_row["POS"], snp_row["P"]] + snp_row.tolist())
    
    # Convert to DataFrames
    columns = ["TRAITID","START", "END", "SNP_POS", "SNP_PVAL"] + tiledb_data.columns.tolist()

    trait_res_df = pd.DataFrame(trait_res, columns=columns).drop(columns=["POS", "P"])
    
    columns = ["TRAITID", "REGION", "SNP_POS", "SNP_PVAL"] + tiledb_data.columns.tolist() + ["S"]
    all_snp_df = pd.DataFrame(all_snp_res, columns=columns).drop(columns=["SNP_POS", "SNP_PVAL"])
    
    return [trait_res_df, all_snp_df]
