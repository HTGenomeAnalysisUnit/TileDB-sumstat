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
from progress.bar import Bar
help_doc = """
Query TileDB database and export data.
"""

@cloup.command("export", no_args_is_help=True, help=help_doc)
@cloup.option_group(
    "Options for querying specific chromosomes, cells, genes or positions in the TileDB",
    cloup.option("--chrom", default = None, type=int, help = "List of chromosomes to filter (e.g. 1,2,3,4)"),
    cloup.option("--cell_file", default = None, type=str, help = "List of cells to interrogate taken from a txt file"),
    cloup.option("--gene_file", default = None, type=str, help = "List of genes taken from a txt file"),
    cloup.option("--snp", default = None, type=str, help = "List of SNPs to interrogate taken from a txt file. Please check README for details on the format of this file"),
    cloup.option("--output_path", default = "out", type=str, help = "Output path with file name where results will be stored")
)

@cloup.option_group(
    "Options for general filters into TileDB",
    cloup.option("--maf", default = 0.01, type=float, help = "The MAF to filter the TILEDB for")
)
@cloup.option_group(
    "Options for Locusbreaker",
    cloup.option("--locusbreaker", is_flag=True, type=bool, help="Option to run locusbreaker"),
    cloup.option("--table", default = None, type=int, help = "Path of the table to provide"),
    cloup.option("--hole-size", default=250000, type=int, help="Minimum pair-base distance between SNPs in different loci (default: 250000)")
)

@click.pass_context
def export(
        ctx,
        uri: str,
        chrom:int,
        cell_file: str,
        gene_file: str,
        snp: str,
        maf: float,
        locusbreaker: bool,
        table: str,
        hole_size: int,
        output_path: str,
        schema: bool
        ):
    
    #Open connection with TileDB
    tiledb_export = tiledb.open(uri, mode="r")
    #Print only the schema of the tiledb
    if schema:
        print(tiledb_export.schema)
        exit()

    #Get list of genes, cell type and positions or create ones
    unique_positions = slice(None)
    if(chrom):
        chrom_list = chrom.split(",")
    else:
        chrom_list = slice(None)
    if(cell_file):
        cell_list = open(cell_file, "r").read().rstrip().split("\n")
    else:
        cell_list = slice(None)
    if(gene_file):
        gene_list = open(gene_file, "r").read().rstrip().split("\n")
    else:
        gene_list = slice(None)

    #Intersect the tiledb with a list of SNPs
    if snp: 
        snp_list = pd.read_table(snp, dtype = {"CHR":str, "POS":np.uint32, "A0":str, "A1":str})
        unique_positions = snp_list['position'].unique().tolist()
        #Open a streaming connection with TileDB
        with tiledb_export as A:
            tiledb_iterator = A.query(
                return_incomplete=True
            ).df[chrom_list, cell_list ,gene_list ,unique_positions]
            #Open a streaming connection with output and run the function
            with open(output_path + ".csv", mode="a") as f:
                for chunk in tiledb_iterator:
                    process_write_chunk(chunk, snp_list, f)
        print(f"Saved filtered summary statistics by SNPs in {output_path}.csv")

    elif locusbreaker:
        print("Starting LocusBreaker")
        tasks = []
        #Defining the Dask functions for delayed
        traits = pd.read_table(table)
        @delayed
        def query_gene(uri, chrom, gene, cell):
            with tiledb.open(uri, mode="r") as tiledb_data:
                return tiledb_data.query(dims=['CHR','CELL','GENE','POS'], attrs=['SNP', 'AF' , 'BETA', 'SE', 'P', 'N']).df[chrom, cell ,gene ,unique_positions]
            
        @delayed
        def delayed_locus_breaker(tiledb_data, pvalue_sig, pvalue_limit, hole_size):
            # Call locus_breaker with the computed tiledb_data
            return locus_breaker(tiledb_data, pvalue_sig=pvalue_sig, pvalue_limit=pvalue_limit, hole_size=hole_size)
        #The computation is divided and run in parallel for each cell and gene separately
        gene_batches = [gene_list[i:i + 100] for i in range(0, len(gene_list), 100)]

        for ind, row in traits.iterrows():
            chrom = row["chrom"]
            cell = row["cell"]
            gene = row["gene"]
            chrom = row["chrom"]
            pvalue_sig = row["pvalue_sig"]
            pvalue_limit = row["pvalue_limit"]
            task = delayed_locus_breaker(query_gene(uri, gene, cell),pvalue_sig=pvalue_sig,pvalue_limit=pvalue_limit,hole_size=hole_size)
            tasks.append(task)

        #The batch size to use which is set to the number of workers if Dask is run
        if ctx.obj["workers"]:
            batch_size = ctx.obj["workers"]
        else:
            batch_size = 1
        computed_results = []
        for i in range(0, len(tasks), batch_size):
                print(f"Batch {i} of {len(tasks)}")
                batch = tasks[i:i+batch_size]
                batch_results = compute(*batch)  # Compute the batch
                for result in batch_results:
                    if not len(result)==0 and not result[0].empty:
                        #print(result)
                        interval = result[0]
                        segments = result[1]
                        interval.to_csv(output_path + "_interval.csv", mode="a", index=False, header = None)
                        segments.to_csv(output_path + "_segment.csv", mode="a", index=False, header = None)
    #If no SNP or locusbreker is run only a filtering is done
    else:
        with tiledb.open(uri, mode="r") as A:
            tiledb_iterator = A.query(
                return_incomplete=True
            ).df[cell_list , gene_list, unique_positions] 
            for chunk in tiledb_iterator:
                chunk.to_csv(output_path + ".csv", mode="a", index=False, header = True)
        print(f"Saved filtered summary statistics in {output_path}")
        exit()
    


