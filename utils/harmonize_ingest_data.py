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
                                usecols = ["variant_id","start_distance","phenotype_id","slope","slope_se","af", "pval_nominal"], 
                                low_memory=False,
                                dtype={"variant_id":str, "start_distance":str, "phenotype_id":str, "slope":np.float32,"slope_se":np.float32, "af":np.float32, "pval_nominal":np.float64}):
        chunk_pl = pl.from_pandas(chunk)
        chunk_pl = chunk_pl.with_columns(
        pl.when(pl.col("start_distance").str.contains("vs"))
        .then(None)
        .otherwise(pl.col("start_distance"))
        .alias("start_distance2")
        )
        vs_count = chunk_pl.filter(pl.col("start_distance2").str.contains("_vs_")).shape[0]
        print(f"Rows with '_vs_': {vs_count}")
        chunk_pl = chunk_pl.with_columns(
        pl.col("start_distance2").cast(pl.Int32).alias("start_distance")
        )
        chrompos_split = chunk_pl.with_columns(
            pl.col("variant_id").str.split_exact("_", 4)
            .struct.rename_fields(["CHROM",'POS','REF','ALT'])
            .alias("fields")
            ).unnest('fields')
    
        # Ensure POS columns are strings
        chrompos_split = chrompos_split.with_columns(pl.col("POS").cast(pl.Utf8))
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
        "cell_type", "gene", "SNP", "start_distance", "position", "allele0", "allele1", "af", "beta", "se", "p-value"
        ])
        dict_type = {"cell_type":"ascii", "position":np.uint32, "SNP":"ascii", "allele0":"ascii", "allele1":"ascii", "start_distance":np.int32, "af":np.float32, "beta":np.float32, "se":np.float32, "p-value":np.float64}
        tiledb.from_pandas(
                        uri=uri,
                        dataframe=chunk_processed.to_pandas(),
                        index_dims=["cell_type", "gene", "position"],
                        column_types=dict_type,
                        mode="append"
                        )

