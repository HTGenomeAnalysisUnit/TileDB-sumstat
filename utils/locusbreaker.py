import pandas as pd
from typing import List
from utils import compute_pheno_variance


def locus_breaker(
    tiledb_data,
    pvalue_sig: float = 5e-8,
    pvalue_limit: float = 5e-6,
    hole_size: int = 250000,
    phenovar: bool=False,
    maf: float = 0.01,
    category: bool = False,

) -> pd.DataFrame:
    """
    Breaking genome in locus
    Returns two dataframes with reiongs and individual SNPs
    :param tiledb_data: TileDBVCF data (default: None)
    :param pvalue_sig: P-value threshold used to create the regions around significant SNPs (default: 5e-8)
    :param pvalue_limit: P-value threshold for loci borders (default: 5e-6)
    :param hole_size: Minimum pair-base distance between SNPs in different loci (default: 250000)
    :param phenovar: If the phenotypic variance need to be computed or not
    :return: DataFrames with the loci information
    """
    # Filter rows based on the p_limit threshold
    tiledb_data = tiledb_data[tiledb_data["P"] < pvalue_limit]
    if phenovar:
        pv = compute_pheno_variance(tiledb_data)
        tiledb_data["S"] = pv
    else:
        tiledb_data["S"] = 1.0


    # If no rows remain after filtering, return an empty DataFrame
    if tiledb_data.empty:
        return []
    # Group by gene first, then calculate regions within each gene
    tiledb_data = tiledb_data[tiledb_data["AF"]>maf]
    tiledb_data = tiledb_data[tiledb_data["AF"]<1-maf]
    if(category=="cis"):
        tiledb_data = tiledb_data[tiledb_data["DIST"]<1000000]
        tiledb_data = tiledb_data[tiledb_data["DIST"]>-1000000]
    elif(category=="trans"):
        tiledb_data = tiledb_data[tiledb_data["DIST"]>1000000]
        tiledb_data = tiledb_data[tiledb_data["DIST"]<-1000000]

    trait_res = []
    trait_res_allsnp = []
    for gene, gene_df in tiledb_data.groupby("GENE"):
        # Find regions where gaps between positions exceed hole_size within each chromosome
        gaps = gene_df["POS"].diff() > hole_size
        group = gaps.cumsum()

        # Group by the identified regions within the chromosome
        for _, group_df in gene_df.groupby(group):
            if group_df["P"].min() < pvalue_sig:
                start_pos = group_df["POS"].min()
                end_pos = group_df["POS"].max()
                best_snp = group_df.loc[group_df["P"].idxmin()]

                region = str(group_df["CHR"].iloc[0]) + ":" + group_df["CELL"].iloc[0]+ ":" + gene + ":" + str(start_pos) + ":" + str(end_pos)


                # Store the interval with the best SNP
                line_res = [gene, start_pos, end_pos, best_snp["POS"], best_snp["P"]] + best_snp.tolist()
                trait_res.append(line_res)

                #Collect all SNPs within the region
                for _, snp_row in group_df.iterrows():
                    snp_res = [region, snp_row["POS"], snp_row["P"]] + snp_row.tolist()
                    trait_res_allsnp.append(snp_res)

    # Convert results to a DataFrame

    columns = ["GENE", "START", "END", "SNP_POS", "SNP_PVAL"] + tiledb_data.columns.tolist()
    trait_res_df = pd.DataFrame(trait_res, columns=columns)
    trait_res_df = trait_res_df.drop(trait_res_df.columns[0], axis=1)
    trait_res_df = trait_res_df.drop(columns=["POS", "P"])
    columns = ["REGION","SNP_POS", "SNP_PVAL"] + tiledb_data.columns.tolist()
    trait_res_allsnp_df = pd.DataFrame(trait_res_allsnp, columns=columns)
    #trait_res_allsnp_df = trait_res_allsnp_df.drop(trait_res_allsnp_df.columns[0], axis=1)
    trait_res_allsnp_df = trait_res_allsnp_df.drop(columns=["SNP_POS", "SNP_PVAL"])

    return [trait_res_df,trait_res_allsnp_df]
