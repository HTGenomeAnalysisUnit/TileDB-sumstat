#Reads the files in the directory and creates a TileDB array from them
import numpy as np
import os
import click 
import cloup
from scqtl.utils.harmonize_ingest import Harmonize
#import dask
import pandas as pd

@cloup.command("ingest", no_args_is_help=True, help="Ingest single cell QTL with TileDB")
@cloup.option_group(
    "Essential parameters",
    cloup.option("--uri-path", default = None, type=str, help = "Where to store the TileDB"),
    cloup.option("--file-path", default = None, type=str, help = "List of the files to ingest"),
    cloup.option("--mapping-file", default = None, type=str, help = "List of the files to ingest"),
    cloup.option("--type-sumstat", default = "scqtl", type=str, help = "Either gwas or scqtl, to specify the type of summary statistics being ingested. This is used to harmonize the data accordingly.")
)
@cloup.option_group(
    "Optional parameters",
    cloup.option("--pvar-file", default = None, type=str, help = "pvar file used to verify the alleles order"),
    cloup.option("--batch-size", default = 1, type=int, help = "The number of files to ingest at once"),
    cloup.option("--chunk-size", default = 50000000, type=int, help = "The number of rows to ingest at once"),
)

def ingest(uri_path: str, mapping_file: str, chunk_size: int, batch_size: int, file_path:str, type_sumstat:str, pvar_file:str = None):
    file_list = pd.read_csv(file_path, sep="\t", header=0, dtype=str)
    #This could be optimized with Dask
    # Create a Harmonize object
    harmonized_object = Harmonize(mapping_file= mapping_file, chunk_size=chunk_size , uri=uri_path, type_sumstat=type_sumstat)
    #CHeck if the tiledb already exists, if not create it
    if not os.path.exists(uri_path):
        print(f"Creating TileDB at {uri_path}")
        harmonized_object.create_tiledb()
    else:
        print(f"TileDB already exists at {uri_path}, skipping creation.")
    harmonized_object.create_mapping()
    cell = None
    gene = None
    trait = None
    N = None
    #for i in range(0, len(file_list), batch_size):
    for record_index, record in file_list.iterrows():
        #batch_files = file_list[i:i + batch_size]
        file = record["FILE"]
        if "N" in file_list.columns:
            N = record["N"]
        if type_sumstat == "scqtl":
            if ["CELL", "GENE"] in file_list.columns:
                cell = record["CELL"]
                gene = record["GENE"]
        if type_sumstat == "gwas":
            if "TRAIT" in file_list.columns:
                trait = record["TRAIT"]
        
        if not os.path.exists(file):
            print(f"File {file} does not exist. Skipping.")
            continue
        print(f"Processing file: {file}")
        # Harmonize the data
        print(f"Harmonizing file: {file}")
        harmonized_object.harmonize(file, trait, cell, gene, N = N)
        # Ingest the data
        print(f"Ingesting data: {file}")
        harmonized_object.ingest_data()
        # Create metadata
        harmonized_object.create_metadata(file_path = file)

        
        #Harmonize(chunk_size, file, uri, pvar_file, type_sumstat="gwas")
        #print(f"Batch {i // batch_size + 1} completed.")   
