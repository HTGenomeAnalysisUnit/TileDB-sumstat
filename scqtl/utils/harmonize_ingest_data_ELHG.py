import numpy as np
import pandas as pd
import scipy.stats as stats
import polars as pl
import tiledb
import glob

def harmonize_ingest_data(chunk_size, file, uri):
    file = file.split(",")
    for chunk in pd.read_table(
        file[0],
        chunksize=int(chunk_size),
        engine="c",
        usecols=[
            "variant_id", "start_distance", "phenotype_id", "slope", "slope_se", "af", "pval_nominal", "ma_count"
        ],
        low_memory=False,
        dtype={
            "variant_id": str,
            "start_distance": str,
            "phenotype_id": str,
            "slope": np.float32,
            "slope_se": np.float32,
            "af": np.float32,
            "pval_nominal": np.float64,
            "ma_count": np.int64
        }
    ):
        chunk_pl = pl.from_pandas(chunk)
        
        chunk_pl = chunk_pl.with_columns(
            pl.when(pl.col("start_distance").str.contains("vs"))
            .then(None)
            .otherwise(pl.col("start_distance"))
            .alias("start_distance2")
        )
        
        chunk_pl = chunk_pl.with_columns(
            pl.col("start_distance2").cast(pl.Int64).alias("DIST"),
            pl.lit(None).alias("RSID")
        )
        
        chrompos_split = chunk_pl.with_columns(
            pl.col("variant_id").str.split_exact("_", 4)
            .struct.rename_fields(["CHR", "POS", "REF", "ALT"])
            .alias("fields")
        ).unnest("fields")
        chrompos_split = chrompos_split.with_columns(pl.col("CHR").str.replace(r"chr", ""))
        
        chrompos_split = chrompos_split.with_columns(
            (pl.col("ma_count") / (pl.col("af") * 2)).alias("N"),
            pl.col("POS").cast(pl.Utf8)
        )
        
        chrom = pl.col("CHR")
        pos = pl.col("POS")
        ref = pl.col("REF")
        alt = pl.col("ALT")
        
        allele0 = pl.when(ref < alt).then(ref).otherwise(alt)
        allele1 = pl.when(ref >= alt).then(ref).otherwise(alt)
        
        sorted_snp = chrom.cast(str) + "_" + pos.cast(str) + "_" + allele0 + "_" + allele1
        new_beta = pl.when(allele0 == ref).then(-pl.col("slope")).otherwise(pl.col("slope"))
        new_af = pl.when(allele0 == ref).then(1 - pl.col("af")).otherwise(pl.col("af"))
        
        chrompos_split = chrompos_split.with_columns([
            pl.lit(file[1]).alias("CELL"),
            pl.col("phenotype_id").alias("GENE"),
            pos.alias("POS"),
            sorted_snp.alias("SNP"),
            new_beta.alias("BETA"),
            pl.col("slope_se").alias("SE"),
            new_af.alias("AF"),
            pl.col("pval_nominal").alias("P")
        ])
        
        chunk_processed = chrompos_split.select([
            "CHR", "CELL", "GENE", "POS", "SNP", "RSID", "DIST", "AF", "BETA", "SE", "P", "N"
        ])
        
        dict_type = {
            "CHR": np.uint16, "CELL": "ascii", "GENE": "ascii", "POS": np.uint32,
            "SNP": "ascii", "RSID": "ascii", "DIST": np.int64, "AF": np.float32,
            "BETA": np.float32, "SE": np.float32, "P": np.float64, "N": np.int64
        }
        
        tiledb.from_pandas(
            uri=uri,
            dataframe=chunk_processed.to_pandas(),
            index_dims=["CHR", "CELL", "GENE", "POS"],
            column_types=dict_type,
            mode="append"
        )
