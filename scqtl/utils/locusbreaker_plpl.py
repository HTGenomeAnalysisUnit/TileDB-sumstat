import polars as pl
from typing import List
from scqtl.utils import compute_pheno_variance

def locusbreaker_plpl(
    tiledb_data,
    pvalue_sig: float = 5e-8,
    pvalue_limit: float = 5e-6,
    hole_size: int = 250000,
    phenovar: bool = False,
    category: bool = False,
    type_sumstat: str = "scqtl",
    maf: float = 0.001,
    locus_max_size: int = 3000000
) -> List[pl.DataFrame]:
    # Convert to Polars DataFrame
    df = pl.from_arrow(tiledb_data)
    
    # Filter by MAF
    df = df.with_columns(
        MAF=pl.when(pl.col("EAF") <= 0.5).then(pl.col("EAF")).otherwise(1 - pl.col("EAF"))
    ).filter(pl.col("MAF") >= maf).drop("MAF")
    
    # Compute phenotypic variance if needed
    if phenovar:
        # Note: compute_pheno_variance might need adjustment for Polars
        df = df.with_columns(S=compute_pheno_variance(df.to_pandas()))
    else:
        df = df.with_columns(S=pl.lit(1.0))
    
    # Filter SNPs for locus definition
    loci_snps = df.filter(pl.col("P") <= pvalue_limit)
    
    if loci_snps.is_empty():
        return []
    
    # Apply cis/trans filtering for scqtl
    if type_sumstat == "scqtl":
        if category == "cis":
            loci_snps = loci_snps.filter(
                (pl.col("DIST") > -1000000) & (pl.col("DIST") < 1000000)
            )
        elif category == "trans":
            loci_snps = loci_snps.filter(
                (pl.col("DIST") < -1000000) | (pl.col("DIST") > 1000000)
            )
    
    # Define grouping keys
    if type_sumstat == "gwas":
        group_keys = ["CHR", "TRAIT"]
    else:
        group_keys = ["CHR", "CELL", "GENE"]
    
    # Process groups
    loci_snps = loci_snps.sort(group_keys + ["POS"])
    
    # Identify groups based on hole_size
    grouped = loci_snps.with_columns(
        pl.col("POS").diff().gt(hole_size).cast(pl.UInt32).fill_null(0).cum_sum().over(group_keys).alias("group_id")
    )
    
    # Aggregate groups
    aggregated = grouped.group_by(group_keys + ["group_id"]).agg(
        [
            pl.col("POS").min().alias("min_pos"),
            pl.col("POS").max().alias("max_pos"),
            pl.col("P").min().alias("min_p"),
            pl.all().sort_by("P").first()
        ]
    )
    
    # Filter significant groups and expand regions
    significant = aggregated.filter(pl.col("min_p") <= pvalue_sig).with_columns(
    start_pos=pl.when(pl.col("min_pos") <= 100000)
            .then(1)
            .otherwise(pl.col("min_pos") - 100000),
            end_pos=pl.col("max_pos") + 100000
    ).filter(pl.col("end_pos") - pl.col("start_pos") < locus_max_size)
    if significant.is_empty():
        return []
    
    hla_start = 28510120
    hla_end = 33480577
    chr17_inv_start = 44849948
    chr17_inv_end = 44899445
    
    # Filter out regions that overlap with excluded regions
    significant = significant.filter(
        ~(
            # Exclude HLA region on chr6
            (pl.col("CHR") == 6) & 
            (
                (pl.col("start_pos") <= hla_end) & (pl.col("end_pos") >= hla_start) |
                (pl.col("start_pos") >= hla_start) & (pl.col("end_pos") <= hla_end) |
                (pl.col("start_pos") <= hla_end) & (pl.col("end_pos") >= hla_end) |
                (pl.col("start_pos") <= hla_start) & (pl.col("end_pos") >= hla_start)
            )
        ) &
        ~(
            # Exclude inversion region on chr17
            (pl.col("CHR") == 17) & 
            (
                (pl.col("start_pos") <= chr17_inv_end) & (pl.col("end_pos") >= chr17_inv_start) |
                (pl.col("start_pos") >= chr17_inv_start) & (pl.col("end_pos") <= chr17_inv_end) |
                (pl.col("start_pos") <= chr17_inv_end) & (pl.col("end_pos") >= chr17_inv_end) |
                (pl.col("start_pos") <= chr17_inv_start) & (pl.col("end_pos") >= chr17_inv_start)
            )
        )
    )

    # Collect all SNPs in significant regionsbj
    regions = significant.select(group_keys + ["start_pos", "end_pos"])
    all_snps = df.join(regions, on=group_keys).filter(
        pl.col("POS").is_between(pl.col("start_pos"), pl.col("end_pos"))
    )
    
    # Format output
    significant = significant.with_columns(
        REGION=pl.format("{}:{}:{}", pl.col("CHR"), pl.col("start_pos"), pl.col("end_pos"))
    )
    
    if type_sumstat == "gwas":
        trait_res_df = significant.select(
            ["TRAIT", "start_pos", "end_pos", "POS", "P"] + [c for c in df.columns if c not in ["POS", "P"]]
        ).rename({"start_pos": "START", "end_pos": "END", "POS": "SNP_POS", "P": "SNP_PVAL"})
        
        all_snp_df = all_snps.with_columns(
            REGION=pl.format("{}:{}:{}", pl.col("CHR"), pl.col("start_pos"), pl.col("end_pos"))
        ).select(
            ["TRAIT", "REGION", "POS", "P"] + [c for c in df.columns if c not in ["POS", "P"]]
        ).rename({"POS": "SNP_POS", "P": "SNP_PVAL"})
        
        trait_res_df = trait_res_df.with_columns(TYPE=pl.lit("gwas"))
        all_snp_df = all_snp_df.with_columns(TYPE=pl.lit("gwas"))
    else:
        trait_res_df = significant.with_columns(
            TRAIT=pl.col("CELL") + ":" + pl.col("GENE")
        ).select(
            ["TRAIT", "start_pos", "end_pos", "POS", "P"] + [c for c in df.columns if c not in ["POS", "P"]]
        ).rename({"start_pos": "START", "end_pos": "END", "POS": "SNP_POS", "P": "SNP_PVAL"})
        
        all_snp_df = all_snps.with_columns(
            TRAIT=pl.col("CELL") + ":" + pl.col("GENE"),
            REGION=pl.format("{}:{}:{}", pl.col("CHR"), pl.col("start_pos"), pl.col("end_pos"))
        ).select(
            ["TRAIT", "REGION", "POS", "P"] + [c for c in df.columns if c not in ["POS", "P"]]
        ).rename({"POS": "SNP_POS", "P": "SNP_PVAL"})
        
        trait_res_df = trait_res_df.with_columns(TYPE=pl.lit("qtl"))
        all_snp_df = all_snp_df.with_columns(TYPE=pl.lit("qtl"))
    
    return [trait_res_df.to_pandas(), all_snp_df.to_pandas()]