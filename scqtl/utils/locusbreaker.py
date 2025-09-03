import pandas as pd
from typing import List
from scqtl.utils import compute_pheno_variance

def locus_breaker(
    tiledb_data,
    pvalue_sig: float = 5e-8,
    pvalue_limit: float = 5e-6,
    hole_size: int = 250000,
    phenovar: bool = False,
    category: bool = False,
    type_sumstat: str = "qtl",
    maf: float = 0.001,
    locus_max_size = 1000000
) -> pd.DataFrame:
    """
    Breaking genome in loci and returning all SNPs within the defined loci boundaries.
    :param tiledb_data: TileDBVCF data
    :param pvalue_sig: P-value threshold for significant SNPs
    :param pvalue_limit: P-value threshold for defining loci
    :param hole_size: Minimum base-pair distance to separate loci
    :param phenovar: Compute phenotypic variance or not
    :param category: cis/trans filter for SNPs based on distance
    :param expansion_size: Number of base pairs to expand loci boundaries
    :param type_sumstat: Type of summary statistics, either "gwas" or "scqtl"
    :return: Two DataFrames, one with loci regions and another with all SNPs in loci
    """
    # Create a copy of the original dataset before filtering    
    # Filter for SNPs below the p-value limit to define loci
    print(f"length before maf {str(maf)} {str(len(tiledb_data))}")
    tiledb_data["MAF"] = tiledb_data["EAF"].where(tiledb_data["EAF"] <= 0.5, 1 - tiledb_data["EAF"])
    tiledb_data = tiledb_data[tiledb_data['MAF'] >= maf]
    tiledb_data = tiledb_data.drop("MAF", axis = 1)
    print(f"length after maf {str(maf)} {str(len(tiledb_data))}")
    if phenovar:
        tiledb_data["S"] = compute_pheno_variance(tiledb_data)
    else:
        tiledb_data["S"] = 1.0
    
    loci_snps = tiledb_data[tiledb_data["P"] < pvalue_limit].copy()
    
    if loci_snps.empty:
        return []
    
    # Apply cis/trans filtering if needed
    if type_sumstat == "qtl":
        if category == "cis":
            loci_snps = loci_snps[(loci_snps["DIST"] > -1000000) & (loci_snps["DIST"] < 1000000)]
        elif category == "trans":
            loci_snps = loci_snps[(loci_snps["DIST"] < -1000000) | (loci_snps["DIST"] > 1000000)]
    
    trait_res = []
    all_snp_res = []
    if type_sumstat == "gwas":
        grouped_loci = loci_snps.groupby(["CHR", "TRAIT"])
    else:
        grouped_loci = loci_snps.groupby(["CHR", "CELL", "GENE"])

    for gene, gene_df in grouped_loci:
        gaps = gene_df["POS"].diff() > hole_size
        group = gaps.cumsum()
        for _, group_df in gene_df.groupby(group):
            if group_df["P"].min() < pvalue_sig:
                
                lower_pos = group_df["POS"].min()
                if int(lower_pos) - 100000 < 0:
                    start_pos = 1
                else:
                    start_pos = lower_pos - 100000
                end_pos = group_df["POS"].max() + 100000
                best_snp = group_df.loc[group_df["P"].idxmin()]
                region = f"{group_df['CHR'].iloc[0]}:{start_pos}:{end_pos}"
                if (end_pos - start_pos) < locus_max_size:
                    trait_res.append([start_pos, end_pos, best_snp["POS"], best_snp["P"]] + best_snp.tolist())
                    # Include all SNPs within the expanded region from the original dataset
                    expanded_snps = tiledb_data[
                        (tiledb_data["CHR"] == group_df["CHR"].iloc[0]) &
                        (tiledb_data["POS"] >= start_pos) &
                        (tiledb_data["POS"] <= end_pos)
                        ]
                    for _, snp_row in expanded_snps.iterrows():
                        all_snp_res.append([region, snp_row["POS"], snp_row["P"]] + snp_row.tolist())

    # Convert to DataFrames
    columns = ["START", "END", "SNP_POS", "SNP_PVAL"] + tiledb_data.columns.tolist()
    if type_sumstat == "gwas":
        trait_res_df = pd.DataFrame(trait_res, columns=columns).drop(columns=["POS", "P"])
        columns = ["REGION", "SNP_POS", "SNP_PVAL"] + tiledb_data.columns.tolist()
        all_snp_df = pd.DataFrame(all_snp_res, columns=columns).drop(columns=["SNP_POS", "SNP_PVAL"])
        trait_res_df = trait_res_df[['TRAIT'] + [col for col in trait_res_df.columns if col != 'TRAIT']]
        all_snp_df = all_snp_df[['TRAIT'] + [col for col in all_snp_df.columns if col != 'TRAIT']]
        trait_res_df[['TYPE']] =  'gwas'
        all_snp_df[['TYPE']] =  'gwas'
    else:
        trait_res_df = pd.DataFrame(trait_res, columns=columns).drop(columns=["POS", "P"])
        columns = ["REGION", "SNP_POS", "SNP_PVAL"] + tiledb_data.columns.tolist() + ["S"]
        all_snp_df = pd.DataFrame(all_snp_res, columns=columns).drop(columns=["SNP_POS", "SNP_PVAL"])
        trait_res_df["TRAIT"] = trait_res_df["CELL"] + ":" + trait_res_df["GENE"]
        all_snp_df["TRAIT"] = all_snp_df["CELL"] + ":" + all_snp_df["GENE"]
        trait_res_df = trait_res_df[['TRAIT'] + [col for col in trait_res_df.columns if col != 'TRAIT']]
        all_snp_df = all_snp_df[['TRAIT'] + [col for col in all_snp_df.columns if col != 'TRAIT']]
        trait_res_df[['TYPE']] =  'qtl'
        all_snp_df[['TYPE']] =  'qtl'

        
    
    return [trait_res_df, all_snp_df]