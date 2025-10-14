#Reads the files in the directory and creates a TileDB array from them
import numpy as np
import os
import click 
import cloup
from scqtl.utils.harmonize_ingest import Harmonize
#import dask
import pandas as pd
import polars as pl

@cloup.command("ingest", no_args_is_help=True, help="Ingest single cell QTL with TileDB")
@cloup.option_group(
    "Essential parameters",
    cloup.option("--uri-path", default = None, type=str, help = "Where to store the TileDB"),
    cloup.option("--file-path", default = None, type=str, help = "List of the files to ingest"),
    cloup.option("--mapping-file", default = None, type=str, help = "Mapping between columns and TileDB types"),
    cloup.option("--type-sumstat", default = "qtl", type=str, help = "Either gwas or qtl, to specify the type of summary statistics being ingested. This is used to harmonize the data accordingly."),
    cloup.option("--type-trait", default = "quant", type=str, help = "In case gwas is chosen this option is either quant or binary.")
)
@cloup.option_group(
    "Optional parameters",
    cloup.option("--pvar-file", default = None, type=str, help = "pvar file used to verify the alleles order"),
    cloup.option("--sep", default = "\t", type=str, help = "pvar file used to verify the alleles order"),
    cloup.option("--chunk-files", is_flag=True, type=bool, default = False, help = "If chunk files in block"),
    cloup.option("--chunk-size", default = 20000000, type=int, help = "The approximate number of rows to ingest at once"),
    cloup.option("--qc", is_flag=True, type=bool, default = False, help = "Harmonize and QC the summary statistics using gwaslab"),
    cloup.option("--only-meta", is_flag=True, type=bool, default = False, help = "Create and ingest metadata")
)

def ingest(uri_path:str, sep:str, type_trait:str, mapping_file:str, chunk_files:bool,  chunk_size:int, file_path:str, type_sumstat:str, pvar_file:str = None, qc:bool = False, only_meta:bool = False):
    file_list = pd.read_csv(file_path, sep=",", header=0, dtype=str)
    # Create a Harmonize object
    harmonized_object = Harmonize(mapping_file= mapping_file, uri=uri_path, type_sumstat=type_sumstat, pvar_file = pvar_file, type_trait = type_trait)
   
    #Check if the tiledb already exists, if not create it
    if not os.path.exists(uri_path):
        print(f"Creating TileDB at {uri_path}")
        harmonized_object.create_tiledb()
    else:
        print(f"TileDB already exists at {uri_path}, skipping creation.")
    harmonized_object.create_mapping()
    cell = None
    gene = None
    trait = None
    n = None
    n_cases = None
    n_controls = None
    buffer_pl = None
    buffer_count = 0
    for record_index, record in file_list.iterrows():
        file = record["FILE"]
        if "N" in file_list.columns:
            n = record["N"]
        if "N_CASES" in file_list.columns:
            n_cases = record["N_CASES"]
        if "N_CONTROLS" in file_list.columns:
            n_controls = record["N_CONTROLS"]       
        if type_sumstat == "qtl":
            if "CELL" in file_list.columns:
                cell = record["CELL"]
            if "GENE" in file_list.columns:
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
        chunk_pl = pl.read_csv(file,separator=sep,low_memory=True ,null_values="NA")
        if chunk_files:
            file_len = len(chunk_pl)
            # Start or extend the buffer
            if buffer_pl is None:
                buffer_pl = chunk_pl
                buffer_count = file_len
            else:
                buffer_pl = pl.concat([buffer_pl, chunk_pl])
                buffer_count += file_len

            # If the buffer reached chunk_size, process it
            if buffer_count >= chunk_size or file in file_list.iloc[-1:]["FILE"].values:
                print(f"Processing buffered chunk of size {buffer_count}")
                harmonized_object.harmonize(
                    sumstat=buffer_pl, trait=trait, cell=cell, gene=gene,
                    n=n, n_cases=n_cases, n_controls=n_controls
                )
                if qc:
                    harmonized_object.qc_sumstat(file_path=file)
                if only_meta:
                    harmonized_object.create_metadata(file_path=file)
                    harmonized_object.ingest_metadata()
                else:
                    harmonized_object.ingest_data(file_path=file)
            
                # Reset buffer
                buffer_pl = None
                buffer_count = 0
            # Continue to next file
            continue
        harmonized_object.harmonize(sumstat = chunk_pl, trait = trait, cell = cell, gene = gene, n = n, n_cases = n_cases, n_controls = n_controls)
        #Performing QC using GWASLAB
        if qc:
            harmonized_object.qc_sumstat(file_path = file)
        # Ingest the data
        if only_meta:
            # Create metadata
            harmonized_object.create_metadata(file_path = file)
            harmonized_object.ingest_metadata()
        else:
            print(f"Ingesting data: {file}")
            harmonized_object.ingest_data(file_path = file)

