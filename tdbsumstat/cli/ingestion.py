#Reads the files in the directory and creates a TileDB array from them
import numpy as np
import os
import click 
import cloup
from tdbsumstat.utils.harmonize_ingest import Harmonize
import pandas as pd
import polars as pl

@cloup.command("ingest", no_args_is_help=True, help="Ingest single cell QTL with TileDB")
@cloup.option_group(
    "Essential parameters",
    cloup.option("--uri-path", default = None, type=str, help = "Where to store the TileDB"),
    cloup.option("--create-tiledb", is_flag=True, type=bool, help = "Create the TileDB before ingestion"),
    cloup.option("--file-path", default = None, type=str, help = "List of the files to ingest"),
    cloup.option("--mapping-file", default = None, type=str, help = "Mapping between columns and TileDB types"),
    cloup.option("--type-sumstat", default = "qtl", type=str, help = "Either gwas or qtl, to specify the type of summary statistics being ingested. This is used to harmonize the data accordingly."),
    cloup.option("--type-trait", default = "quant", type=str, help = "In case gwas is chosen this option is either quant or binary.")
)
@cloup.option_group(
    "Optional parameters",
    cloup.option("--pvar-file", default = None, type=str, help = "pvar file used to verify the alleles order"),
    cloup.option("--sep", default = "\t", type=str, help = "pvar file used to verify the alleles order"),
    cloup.option("--qc", is_flag=True, type=bool, default = False, help = "Harmonize and QC the summary statistics using gwaslab"),
    cloup.option("--mac", default=None, type=int, help = "Minor allele count filter during ingestion"),
    cloup.option("--permuted", is_flag=True, type=bool, default = False, help = "Compute the SE from the permuted pvalue"),
    cloup.option("--only-meta", is_flag=True, type=bool, default = False, help = "Create and ingest metadata")
)

def ingest(uri_path:str, 
           sep:str, 
           type_trait:str, 
           mapping_file:str, 
           file_path:str, 
           create_tiledb:bool, 
           type_sumstat:str, 
           pvar_file:str = None, 
           qc:bool = False, 
           only_meta:bool = False,
           mac:int = None,
           permuted:bool = False):
    # Create a Harmonize object
    print("Starting ingestion")
    harmonized_object = Harmonize(mapping_file= mapping_file, uri=uri_path, type_sumstat=type_sumstat, pvar_file = pvar_file, type_trait = type_trait, mac = mac)
   
    #Check if the tiledb already exists, if not create it
    if create_tiledb:
        print(f"Creating TileDB at {uri_path}")
        harmonized_object.create_tiledb()
        return
    
    if only_meta:
        print("Merging individual metadata files into final metadata")
        harmonized_object.merge_metadata_files()
        return

    file_list = pd.read_csv(file_path, sep=",", header=0, dtype=str)
    harmonized_object.create_mapping()

    for record_index, record in file_list.iterrows():
        file = record["FILE"]
        n = record.get("N")
        n_cases = record.get("N_CASES")
        n_controls = record.get("N_CONTROLS")
        cell = record.get("CELL")
        gene = record.get("GENE")
        pheno_var = record.get("PHENO_VAR", 1)
        trait = record.get("TRAIT")
        
        if not os.path.exists(file):
            print(f"File {file} does not exist. Skipping.")
            continue
        print(f"Processing file: {file}")
        # Harmonize the data
        print(f"Harmonizing file: {file}")
        chunk_pl = pl.read_csv(file,separator=sep,low_memory=True ,null_values="NA")
        harmonized_object.harmonize(sumstat = chunk_pl, trait = trait, cell = cell, 
                                    gene = gene, pheno_var = pheno_var, n = n, n_cases = n_cases, 
                                    n_controls = n_controls, mac = mac, permuted = permuted)
        #Performing QC using GWASLAB
        if qc:
            harmonized_object.qc_sumstat(file_path = file)
        # Ingest the data
        else:
            print(f"Ingesting data: {file}")
            harmonized_object.ingest_data(file_path = file)
            harmonized_object.create_metadata(file_path = file)
