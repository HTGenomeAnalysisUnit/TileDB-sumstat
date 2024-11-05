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
    cloup.option("--cell_type", default = None, type=str, help = "list of cells to interrogate from a txt file"),
    cloup.option("--genes", default = None, type=str, help = "list of genes from a txt file"),
    cloup.option("--positions", default = None, type=int, help = "start position of the region to query"),
    cloup.option("--SNP_list", default = None, type=str, help = "list of regions to interrogate in the format start-end taken from a txt file"),
    cloup.option("--pvalue-sig", default = None, type=str, help = "output folder where queries will be stored"),
    cloup.option("--out", default = "out", type=str, help = "output folder where queries will be stored"),
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
        positions: str,
        SNP_list: str,
        out: str,
        locusbreaker: bool,
        pvalue_sig : float,
        pvalue_limit : float,
        hole_size : int
        ):
    
    tiledb_s = tiledb.open(uri, mode="r")
    cell_list = open(cells, "r").read().rstrip().split("\n")
    gene_list = open(genes, "r").read().rstrip().split("\n")
    
    if(len(gene_list)>100):
        print("please give a number of genes to query not over 100")
        return
    
    if(SNP_list):
        SNP_list = pd.read_csv(SNP_list)
        position_list = SNP_list[["position"]].to_list()
        subset_SNPs = tiledb_s.query(dims=['cell_type','position','gene'], attrs=['SNP', 'beta', 'p-value']).df[cell_list, gene_list ,position_list]
        subset_SNPs = subset_SNPs.merge(SNP_list, on = "position")

    if positions:
        position_range= positions.split("-")
        start = int(position_range[0])
        end = int(position_range[1])
        if((end - start) > 20000000):
            print("region to query is too big, please provide a smaller region")
            return
        subset_SNPs = tiledb_s.query(dims=['cell_type','position','gene'], attrs=['SNP', 'beta', 'p-value']).df[cell_list, gene_list ,start:end]

    if locusbreaker:
        tasks = []
        @delayed
        def query_gene(tiledb_data, gene, cell):
            return tiledb_s.query(return_arrow = True, dims=['position'], attrs=['SNP', 'beta','p-value']).df[cell, gene, :].to_pandas()
        for cell in cell_list:
            for gene in gene_list:
                task = locus_breaker(query_gene(tiledb_s,gene,cell), out = f"{out}_{cell}_{gene}.parquet")  # Create a delayed task for each gene
                tasks.append(task)
        
        batch_size = 50  # Example batch size
        computed_results = []
        for i in range(0, len(tasks), batch_size):
            batch = tasks[i:i+batch_size]
            batch_results = compute(*batch)  # Compute the batch
