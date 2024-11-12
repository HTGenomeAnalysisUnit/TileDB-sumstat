import tiledb
import click
import cloup
import pandas as pd
import pyarrow.parquet
from dask import delayed, compute
from utils.process_write_chunk import process_write_chunk
from utils.locusbreaker import locus_breaker
import numpy as np
from pyarrow import csv
help_doc = """
Query TileDB database and export data.
"""

@cloup.command("export", no_args_is_help=True, help=help_doc)
@cloup.option_group(
    "Options for querying the TileDB",
    cloup.option("--uri", default = None, type=str, help = "Where to data to be created or queried is stored"),
    cloup.option("--schema", is_flag = True, type=bool, help = "Print the schema of a tiledb"),
    cloup.option("--cell_types", default = None, type=str, help = "List of cells to interrogate taken from a txt file"),
    cloup.option("--genes", default = None, type=str, help = "List of genes taken from a txt file"),
    cloup.option("--snp", default = None, type=str, help = "List of SNPs to interrogate taken from a txt file. Please check README for details on the format of this file"),
    cloup.option("--genes", default = None, type=str, help = "List of genes taken from a txt file"),
    cloup.option("--output_path", default = "out", type=str, help = "Output path with file name where results will be stored"),
)
@cloup.option_group(
    "Options for Locusbreaker",
    cloup.option("--locusbreaker", is_flag=True, type=bool, help="Option to run locusbreaker"),
    cloup.option("--pvalue-sig", default=5e-8, type=float, help="P-value threshold to use for filtering the data"),
    cloup.option("--pvalue-limit", default=5e-6, type=float, help="P-value threshold for loci borders"),
    cloup.option("--hole-size", default=250000, type=int, help="Minimum pair-base distance between SNPs in different loci (default: 250000)")
)

@click.pass_context
def export(
        ctx,
        uri: str,
        cell_types: str,
        genes: str,
        snp: str,
        locusbreaker: bool,
        pvalue_sig: float,
        pvalue_limit:float,
        hole_size: int,
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
    if snp: 
        snp_list = pd.read_table(snp, dtype = {"chr":str, "position":np.uint32, "A0":str, "A1":str})
        unique_positions = snp_list['position'].unique().tolist()
        with tiledb.open(uri, mode="r") as A:
            tiledb_iterator = A.query(
                return_incomplete=True
            ).df[cell_list ,gene_list ,unique_positions]  # Replace with appropriate filters if necessary

            with open(output_path + ".csv", mode="a") as f:
                for chunk in tiledb_iterator:
                    # Convert the chunk to Polars format for processing
                    process_write_chunk(chunk, snp_list, f)
        print(f"Saved filtered summary statistics by SNPs in {output_path}.csv")
    
    elif locusbreaker:
        tasks = []
        @delayed
        def query_gene(tiledb_data, gene, cell):
            return tiledb_s.query(dims=['cell_type','gene','position'], attrs=['SNP' ,'allele0', 'allele1' ,'af' , 'beta', 'se', 'p-value']).df[cell, gene, unique_positions]
        @delayed
        def delayed_locus_breaker(tiledb_data, pvalue_sig, pvalue_limit, hole_size):    
            # Call locus_breaker with the computed tiledb_data
            return locus_breaker(tiledb_data, pvalue_sig=pvalue_sig, pvalue_limit=pvalue_limit, hole_size=hole_size)

        tiledb_s = tiledb.open(uri, mode="r")
        gene_arrow = tiledb_s.query(return_arrow = True, dims=['gene']).df[cell_list, gene_list, unique_positions]
        gene_array = gene_arrow['gene']
        unique_genes = list(set(gene_array.to_pylist()))
        for cell in cell_list:
            for gene in unique_genes:
                task = delayed_locus_breaker(query_gene(tiledb_s, gene, cell),pvalue_sig=pvalue_sig,pvalue_limit=pvalue_limit,hole_size=hole_size)
                tasks.append(task)
        
        batch_size = 50  # Example batch size
        computed_results = []
        for i in range(0, len(tasks), batch_size):
            batch = tasks[i:i+batch_size]
            batch_results = compute(*batch)  # Compute the batch
            
            for result in batch_results:
                interval = result[0]
                segments = result[1]
                interval.to_csv(output_path + "_interval.csv", mode="a", index=False)
                segments.to_csv(output_path + "_segment.csv", mode="a", index=False)
                #batch_results.to_csv(output_path, mode = "a", index = False)

    else:
        with tiledb.open(uri, mode="r") as A:
            tiledb_iterator = A.query(
                return_incomplete=True
            ).df[cell_list , gene_list, unique_positions] 
            for chunk in tiledb_iterator:
                chunk.to_csv(output_path + ".csv", mode="a", index=False, header = True)
        print(f"Saved filtered summary statistics in {output_path}")
        exit()
    


