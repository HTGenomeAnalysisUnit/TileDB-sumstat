import tiledb
import click
import cloup
import pandas as pd
import pyarrow.parquet
from dask import delayed, compute
from utils.process_write_chunk import process_write_chunk
import numpy as np
from pyarrow import csv
help_doc = """
Qeury TileDB database and export data.
"""

@cloup.command("export", no_args_is_help=True, help=help_doc)
@cloup.option_group(
    "Options for querying the TileDB",
    cloup.option("--uri", default = None, type=str, help = "Where to data to be created or queried is stored"),
    cloup.option("--schema", is_flag = True, type=bool, help = "Print the schema of a tiledb"),
    cloup.option("--cell_types", default = None, type=str, help = "List of cells to interrogate taken from a txt file"),
    cloup.option("--genes", default = None, type=str, help = "List of genes taken from a txt file"),
    cloup.option("--snp", default = None, type=str, help = "List of SNPs to interrogate taken from a txt file. Please check README for details on the format of this file"),
    cloup.option("--output_path", default = "out.csv", type=str, help = "Output path with file name where results will be stored"),
)

@click.pass_context
def export(
        ctx,
        uri: str,
        cell_types: str,
        genes: str,
        snp: str,
        output_path: str,
        schema: bool
        ):
    
    #Define an empty slice in case cell_types, genes or positions are not indicated
    cell_list = slice(None)
    gene_list = slice(None)
    unique_positions = slice(None)

    if cell_types:
        cell_list = open(cell_types, "r").read().rstrip().split("\n")
    if(genes):
        gene_list = open(genes, "r").read().rstrip().split("\n")

    tiledb_export = tiledb.open(uri, mode="r")

    #Print only the schema of the tiledb
    if schema:
        print(tiledb_export.schema)
        exit()

    #Intersect the tiledb with a list of SNPs
    if(snp): 
        snp_list = pd.read_table(snp, dtype = {"chr":str, "position":np.uint32, "A0":str, "A1":str})
        unique_positions = snp_list['position'].unique().tolist()
        with tiledb.open(uri, mode="r") as A:
            tiledb_iterator = A.query(
                return_incomplete=True
            ).df[cell_list ,gene_list ,unique_positions]  # Replace with appropriate filters if necessary

            with open(output_path, mode="a") as f:
                for chunk in tiledb_iterator:
                    # Convert the chunk to Polars format for processing
                    process_write_chunk(chunk, snp_list, f)
        print(f"Saved filtered summary statistics by SNPs in {output_path}")

    else:
        with tiledb.open(uri, mode="r") as A:
            tiledb_iterator = A.query(
                return_incomplete=True
            ).df[cell_list , gene_list, unique_positions] 
            for chunk in tiledb_iterator:
                chunk.to_csv(output_path, mode="a", index=False)
        print(f"Saved filtered summary statistics in {output_path}")
        exit()
    


