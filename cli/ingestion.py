#Reads the files in the directory and creates a TileDB array from them
import tiledb
import numpy as np
import os
import click 
import cloup
import pandas as pd
from utils.create_tiledb_schema import create_tiledb_schema
from utils.harmonize_data import harmonize_data

@cloup.command("ingestion", no_args_is_help=True, help="Ingest single cell QTL with TileDB")
@cloup.option_group(
    "Essential parameters",
    cloup.option("--uri", default = None, type=str, help = "Where to store the TileDB"),
    cloup.option("--input", default = None, type=str, help = "A tsv file to ingest"),
    cloup.option("--cell_type", default = None, type=str, help = "The celltype to ingest")
)
@cloup.option_group(
    "Optional parameters",
    cloup.option("--chunk_size", default = 50000000, type=str, help = "The number of rows to ingest at once"),
    
)
@click.pass_context
def ingestion(ctx, uri: str, input: str, chunk_size: int, cell_type: str):
    ctx.obj = {"uri": uri, "input": input, "chunksize": chunk_size}
    if not os.path.exists(ctx.obj["uri"]):
        create_tiledb_schema(ctx.obj["uri"])
    for chunk in pd.read_table(input, 
                                chunksize=int(chunk_size), 
                                compression="gzip", 
                                engine = "c", 
                                usecols = ["variant_id","phenotype_id","slope","slope_se","af", "pval_nominal"], 
                                low_memory=False,
                                dtype={"variant_id":str, "phenotype_id":str, "slope":np.float32,"slope_se":np.float32, "af":np.float32, "pval_nominal":np.float64}):
                    chunk_harmonized = harmonize_data(chunk, cell_type)
                    dict_type = {"cell_type":"ascii", "position":np.uint32, "SNP":"ascii", "allele0":"ascii", "allele1":"ascii", "af":np.float32, "beta":np.float32, "se":np.float32, "p-value":np.float64}
                    tiledb.from_pandas(
                    uri=uri,
                    dataframe=chunk_harmonized,
                    index_dims=["cell_type", "gene", "position"],
                    column_types=dict_type,
                    mode="append"
                    )
                                
    # Ingest data into TileDB
   
