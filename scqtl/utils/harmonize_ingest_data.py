import numpy as np
import pandas as pd
import scipy.stats as stats
import polars as pl
import tiledb

def harmonize_ingest_data(chunk_size, file, uri):
    file = file.split(",")
    for chunk in pd.read_csv(file[0], 
                                chunksize=int(chunk_size), 
                                engine = "c", 
                                usecols = ["chr", "pos", "ref", "variant_id", "alt","start_distance","phenotype_id","slope","slope_se","af", "pval_nominal", "N"], 
                                low_memory=False,
                                dtype={"chr":np.uint16, "pos": np.uint32, "variant_id":str, "start_distance":str, "phenotype_id":str, "slope":np.float32,"slope_se":np.float32, "af":np.float32, "pval_nominal":np.float64, "N":np.int64}):
        chunk_pl = pl.from_pandas(chunk)
        chunk_pl = chunk_pl.with_columns(
        pl.when(pl.col("start_distance").str.contains("vs"))
        .then(None)
        .otherwise(pl.col("start_distance"))
        .alias("start_distance2")
        )
        chunk_pl = chunk_pl.with_columns(
        pl.col("start_distance2").cast(pl.Int64).alias("DIST")
        )
        
        chrom = pl.col("chr")
        pos = pl.col("pos")
        ref = pl.col("ref")
        alt = pl.col("alt")

        # Create computed columns using expressions
        allele0 = pl.when(ref < alt).then(ref).otherwise(alt)
        allele1 = pl.when(ref >= alt).then(ref).otherwise(alt)

        # Create sorted SNP and new beta values using conditional expressions
        sorted_snp = chrom.cast(str) + "_" + pos.cast(str) + "_" + allele0 + "_" + allele1
        new_beta = pl.when(allele0 == ref).then(-pl.col("slope")).otherwise(pl.col("slope"))
        new_af = pl.when(allele0 == ref).then(1 - pl.col("af")).otherwise(pl.col("af"))

        # Add new columns to the DataFrame with expressions
        chrompos_split = chunk_pl.with_columns([
            pl.col("chr").alias("CHR"),
            pl.lit(file[1]).alias("CELL"),  # Assuming `cell_type` is defined somewhere
            pl.col("phenotype_id").alias("GENE"),
            pos.alias("POS"),
            sorted_snp.alias("SNP"),
            pl.col("variant_id").alias("RSID"),
            new_beta.alias("BETA"),
            pl.col("slope_se").alias("SE"),
            new_af.alias("AF"),
            pl.col("pval_nominal").alias("P")
        ])       

        # Select necessary columns for the output
        chunk_processed = chrompos_split.select([
        "CHR", "CELL", "GENE", "POS", "SNP", "RSID", "DIST", "AF", "BETA", "SE", "P","N"
        ])
        dict_type = {"CHR":np.uint16, "CELL":"ascii", "GENE": "ascii", "POS":np.uint32, "SNP":"ascii", "RSID":"ascii", "DIST":np.int64, "AF":np.float32, "BETA":np.float32, "SE":np.float32, "P":np.float64, "N":np.int64}
        tiledb.from_pandas(
                        uri=uri,
                        dataframe=chunk_processed.to_pandas(),
                        index_dims=["CHR", "CELL", "GENE", "POS"],
                        column_types=dict_type,
                        mode="append"
                        )
        with tiledb.open(uri, mode="r+") as A:
            if "GENE_CELLTYPE" not in A.meta:
                A.meta["GENE_CELLTYPE"] = {file[1]: chunk_pl["phenotype_id"].unique().to_list()}
            else:
                if file[1] not in A.meta["GENE_CELLTYPE"]:
                    A.meta["GENE_CELLTYPE"][file[1]] = chunk_pl["phenotype_id"].unique().to_list()
                else:
                    # Append the new cell type to the existing list
                    A.meta["GENE_CELLTYPE"][file[1]].extend(chunk_pl["phenotype_id"].unique().to_list())

            #The line on top is good but i need to append in case the cell already exists

