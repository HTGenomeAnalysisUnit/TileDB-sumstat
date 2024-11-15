import numpy as np
import pandas as pd
import scipy.stats as stats
import polars as pl
import tiledb

def harmonize_ingest_data(chunk_size, file, uri):
    file = file.split(",")
    for chunk in pd.read_table(file[0], 
                                chunksize=int(chunk_size), 
                                compression="gzip", 
                                engine = "c", 
                                usecols = ["variant_id","phenotype_id","slope","slope_se","af", "pval_nominal"], 
                                low_memory=False,
                                dtype={"variant_id":str, "phenotype_id":str, "slope":np.float32,"slope_se":np.float32, "af":np.float32, "pval_nominal":np.float64}):
        #mapping_SNP = pd.read_table(pvar_map, dtype = {"POS":np.int64})
        chunk_pl = pl.from_pandas(chunk)
        chrompos_split = chunk_pl.with_columns(
            pl.col("variant_id").str.split_exact("_", 4)
            .struct.rename_fields(["CHROM",'POS','REF','ALT'])
            .alias("fields")
            ).unnest('fields')
    

        # Ensure POS columns are strings
        chrompos_split = chrompos_split.with_columns(pl.col("POS").cast(pl.Utf8))
        #mapping_SNP = pl.from_pandas(mapping_SNP)
        #mapping_SNP = mapping_SNP.with_columns(pl.col("POS").cast(pl.Utf8))

        # Perform inner join
        #chrompos_split_inner = chrompos_split.join(mapping_SNP, on=["CHROM", "POS", "REF", "ALT"], how="inner")

        # Perform anti-join directly
        #chrompos_split_antijoin = chrompos_split.join(mapping_SNP, on=["CHROM", "POS", "REF", "ALT"], how="anti")

        # Reorder columns for anti-join
        #chrompos_split_outer_fix = chrompos_split_antijoin.with_columns([
        #pl.col("REF").alias("TEMP"),    # Rename "REF" temporarily to "TEMP"
        #pl.col("ALT").alias("REF")      # Rename "ALT" to "REF"
        #]).drop("ALT")                      # Drop the original "ALT" column before renaming "TEMP"

        # Rename "TEMP" to "ALT" after dropping the duplicate
        #chrompos_split_outer_fix = chrompos_split_outer_fix.rename({"TEMP": "ALT"})
    
        # Concatenate the results
        #chrompos_split = pl.concat([chrompos_split_outer_fix, chrompos_split_inner])

        #chrompos_split = chrompos_split_merged.to_pandas()
        chrom = pl.col("CHROM")
        pos = pl.col("POS")
        ref = pl.col("REF")
        alt = pl.col("ALT")

        # Create computed columns using expressions
        allele0 = pl.when(ref < alt).then(ref).otherwise(alt)
        allele1 = pl.when(ref >= alt).then(ref).otherwise(alt)

        # Create sorted SNP and new beta values using conditional expressions
        sorted_snp = chrom + "_" + pos.cast(str) + "_" + allele0 + "_" + allele1
        new_beta = pl.when(allele0 == ref).then(-pl.col("slope")).otherwise(pl.col("slope"))
        new_af = pl.when(allele0 == ref).then(1 - pl.col("af")).otherwise(pl.col("af"))

        # Add new columns to the DataFrame with expressions
        chrompos_split = chrompos_split.with_columns([
            sorted_snp.alias("SNP"),
            pl.col("phenotype_id").alias("gene"),
            pos.alias("position"),
            ref.alias("allele0"),
            alt.alias("allele1"),
            new_beta.alias("beta"),
            pl.col("slope_se").alias("se"),
            pl.lit(file[1]).alias("cell_type"),  # Assuming `cell_type` is defined somewhere
            new_af.alias("af"),
            pl.col("pval_nominal").alias("p-value")
        ])

        # Select necessary columns for the output
        chunk_processed = chrompos_split.select([
        "cell_type", "gene", "SNP", "position", "allele0", "allele1", "af", "beta", "se", "p-value"
        ])
        dict_type = {"cell_type":"ascii", "position":np.uint32, "SNP":"ascii", "allele0":"ascii", "allele1":"ascii", "af":np.float32, "beta":np.float32, "se":np.float32, "p-value":np.float64}
        tiledb.from_pandas(
                        uri=uri,
                        dataframe=chunk_processed.to_pandas(),
                        index_dims=["cell_type", "gene", "position"],
                        column_types=dict_type,
                        mode="append"
                        )

