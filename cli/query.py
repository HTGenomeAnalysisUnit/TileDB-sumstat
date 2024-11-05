import tiledb
import click
import cloup
import pandas as pd
import pyarrow.parquet
from methods.locusbreaker import locus_breaker
from dask import delayed, compute

help_doc = """
Qeury TileDB data by cell types, genes and regions.
"""

@cloup.command("query", no_args_is_help=True, help=help_doc)
@cloup.option_group(
    "Options for querying the TileDB",
    cloup.option("--uri", default = None, type=str, help = "Where to data to be created or queried is stored"),
    cloup.option("--cells", default = None, type=str, help = "list of cells to interrogate from a txt file"),
    cloup.option("--genes", default = None, type=str, help = "list of genes from a txt file"),
    cloup.option("--start", default = None, type=int, help = "start position of the region to query"),
    cloup.option("--end", default = None, type=int, help = "end position of the region to query"),
    cloup.option("--region_list", default = None, type=str, help = "list of regions to interrogate in the format start-end taken from a txt file"),
    cloup.option("--out", default = "out", type=str, help = "output folder where queries will be stored")
)
@cloup.option_group(
    "Options for Locusbreaker",
    cloup.option("--locusbreaker", default=False, is_flag=True, type=bool, help="Option to run locusbreaker"),
    cloup.option("--pvalue-sig", default=5e-5, type=float, help="P-value threshold to use for filtering the data"),
    cloup.option("--pvalue-limit", default=5e-5, type=float, help="P-value threshold for loci borders"),
    cloup.option("--hole-size", default=250000, type=int, help="Minimum pair-base distance between SNPs in different loci (default: 250000)")
)
@click.pass_context
def query(
        ctx,
        uri: str,
        cells: str,
        genes: str,
        start: int,
        end: int,
        region_list: str,
        out: str,
        locusbreaker: bool,
        pvalue_sig : float,
        pvalue_limit : float,
        hole_size : int
        ):

    
    l_cells = open(cells, "r").read().rstrip().split("\n")
    l_genes = open(genes, "r").read().rstrip().split("\n")
    if(len(l_genes)>100):
        print("please give a number of genes to query not over 100")
        return
    if locusbreaker:
        tiledb_s = tiledb.open(uri, mode="r")
        tasks = []
        @delayed
        def query_gene(tiledb_data, gene, cell):
            return tiledb_s.query(return_arrow = True, dims=['position'], attrs=['SNP', 'beta','p-value']).df[cell, gene, :].to_pandas()
        for cell in l_cells:
            for gene in l_genes:
                task = locus_breaker(query_gene(tiledb_s,gene,cell), out = f"{out}_{cell}_{gene}.parquet")  # Create a delayed task for each gene
                tasks.append(task)
        
        batch_size = 50  # Example batch size
        computed_results = []
        for i in range(0, len(tasks), batch_size):
            batch = tasks[i:i+batch_size]
            batch_results = compute(*batch)  # Compute the batch

    else:
        if((end - start) > 20000000):
            print("region to query is too big, please provide a smaller region")
            return
        else:
            tiledb_s = tiledb.open(uri, mode="r")
            tiledb_q = tiledb_s.query(return_arrow=True, dims = ["cell_type", "gene"], attrs = ["SNP", "beta", "p-value"]).df[l_cells, l_genes , start:end]
            pyarrow.parquet.write_table(tiledb_q, f"{out}.parquet")

