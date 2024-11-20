#Reads the files in the directory and creates a TileDB array from them
import numpy as np
import os
import click 
import cloup
from utils.create_tiledb_schema import create_tiledb_schema
from utils.harmonize_ingest_data import harmonize_ingest_data
import dask

@cloup.command("ingestion", no_args_is_help=True, help="Ingest single cell QTL with TileDB")
@cloup.option_group(
    "Essential parameters",
    cloup.option("--uri", default = None, type=str, help = "Where to store the TileDB"),
    cloup.option("--list_files", default = None, type=str, help = "List of the files to ingest")
)
@cloup.option_group(
    "Optional parameters",
    cloup.option("--batch_size", default = 50000000, type=int, help = "The number of rows to ingest at once")
)
@click.pass_context
def ingestion(ctx, uri: str, chunk_size: int, batch_size: int, list_files:str):
    #ctx.obj = {"uri": uri, "input": input, "chunksize": chunk_size}
    if not os.path.exists(uri):
        create_tiledb_schema(uri)
    file_list = open(list_files, "r").read().splitlines()
    for i in range(0, len(file_list), batch_size):
        batch_files = file_list[i:i + batch_size]
        tasks = [dask.delayed(harmonize_ingest_data)(chunk_size, file, uri) for file in batch_files]
        dask.compute(*tasks)
        print(f"Batch {i // batch_size + 1} completed.")   
