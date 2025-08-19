import tiledb
import click
import cloup
import pandas as pd
from dask import delayed, compute
#from scqtl.utils.process_write_chunk import process_write_chunk
from scqtl.utils.locusbreaker import locus_breaker
import numpy as np
import os

help_doc = """
Query TileDB database and export data.
"""

@cloup.command("export", no_args_is_help=True, help=help_doc)
@cloup.option_group(
    "Options for querying specific chromosomes, cells, genes or positions in the TileDB",
    cloup.option("--uri-path", default = None, type=str, help = "path of TileDB"),
    cloup.option("--chrom", default = None, type=int, help = "chromosome to filter (e.g. 1,2,3,4)"),
    cloup.option("--trait", default = None, type=str, help = "Cell to interrogate"),
    cloup.option("--cell", default = None, type=str, help = "Cell to interrogate"),
    cloup.option("--gene", default = None, type=str, help = "Genes to interrogate"),
    cloup.option("--table_regions", default = None, type=str, help = "Regions to interrogate from a table"),
    cloup.option("--attr", default = "P,SNPID,EAF,BETA,SE,N", type=str, help = "Attributes to output"),
    cloup.option("--snp", default = None, type=str, help = "List of SNPs to interrogate taken from a txt file. Please check README for details on the format of this file")
)

@cloup.option_group(
    "Options for getting info from TileDB",
    cloup.option("--schema", is_flag=True, type=bool, default = False, help="Option to print on screen the schema of the TileDB used"),
)
@cloup.option_group(
    "Options for Locusbreaker",
    cloup.option("--locusbreaker", is_flag=True, type=bool, default = False, help="Option to run locusbreaker"),
    cloup.option("--pvalue_sig", default = 5e-8, type=float, help = "P-value threshold used to create the regions around significant SNPs (default: )"),
    cloup.option("--pvalue_limit", default = 1e-5, type=float, help = "P-value threshold for loci borders"),
    cloup.option("--hole", default = 250000, type=int, help = "Minimum pair-base distance between SNPs in different loci (default: 250000)"),
    cloup.option("--phenovar", is_flag = True, type=bool, default = False, help = "Compute the phenotypic variance"),
    cloup.option("--maf", default = 0.01, type=float, help = "The MAF to filter the TILEDB before locusbreaker is run"),
    cloup.option("--category",default = "cis",type=str,  help = "If locusbreaker run on cis or trans QLTs"),
    cloup.option("--table", default = None, type=str, help = "Path of the table to provide"),
    cloup.option("--type-sumstat", default = None, type=str, help = "Path of the table to provide"),
)
@cloup.option_group(
    "Options for output",
    cloup.option("--out_lb", default = "out", type=str, help = "Output path with file name where results will be stored"),
    cloup.option("--out_rg", default = "out", type=str, help = "Output path with file name where results will be stored")
)

@click.pass_context
def export(
        ctx,
        uri_path: str,
        schema: bool,
        chrom:int,
        cell: str,
        gene: str,
        type_sumstat: str,
        trait: str,
        table_regions: str,
        attr: str,
        snp: str,
        locusbreaker: bool,
        maf: float,
        category: str,
        phenovar: bool,
        table: str,
        pvalue_sig: float,
        pvalue_limit: float,
        hole: int,
        out_lb: str,
        out_rg: str
        ):
    
    #Open connection with TileDB
    tiledb_export = tiledb.open(uri_path, mode="r")
    client = ctx.obj.get("dask_cluster", None)
    #Print only the schema of the tiledb

    #Get list of genes, cell type and positions or create ones
    unique_positions = slice(None)
    if not chrom:
        chrom = slice(None)
    if not cell:
        cell = slice(None)
    if not gene:
        gene = slice(None)
    if not trait:
        trait = slice(None)

    if schema:
        print(tiledb_export.schema)
        if client:
            print("Shutting down Dask cluster...")
            client.close()
        exit()

    #Intersect the tiledb with a list of SNPs
    if snp: 
        snp_list = pd.read_table(snp, dtype = {"CHR":str, "POS":np.uint32, "A0":str, "A1":str})
        unique_positions = snp_list['position'].unique().tolist()
        #Open a streaming connection with TileDB
        for chrom in chrom_list:
            if type_sumstat == "gwas":
                chrom_list = snp_list['CHR'].unique().tolist()
                with tiledb_export as A:
                    tiledb_iterator = A.query(
                        return_incomplete=True,
                        attrs=attr.split(",")
                    ).df[chrom, trait ,unique_positions]
            else:
                chrom_list = snp_list['CHR'].unique().tolist()
                with tiledb_export as A:
                    tiledb_iterator = A.query(
                        return_incomplete=True,
                        attrs=attr.split(",")
                    ).df[chrom, cell , gene, unique_positions]
            #Open a streaming connection with output and run the function
            with open(out_rg + ".csv", mode="a") as f:
                for chunk in tiledb_iterator:
                    chunk.to_csv(out_rg + ".csv", mode="a", index=False, header = True)
        print(f"Saved filtered summary statistics by SNPs in {out_rg}.csv")
        client.close()
    elif table_regions:
        pd_region = pd.read_csv(table_regions)
        counter_nonempty_region = 0
        for ind, row  in pd_region.iterrows():
            if type_sumstat == "gwas":
                trait = row["TRAIT"]
                region = tiledb_export.query(dims = ["CHR","POS", "TRAIT"], attrs = attr.split(",")).df[int(row["CHR"]),trait,int(row["START"]):int(row["END"])]
            
            else:
                cell,gene = row["TRAIT"].split(":")
                region = tiledb_export.query(dims = ["CHR","POS","CELL","GENE"], attrs = attr.split(",")).df[int(row["CHR"]),cell,gene,int(row["START"]):int(row["END"])]
            if len(region)>0:
                if counter_nonempty_region==0:
                    region.to_csv(out_rg, mode='a', index = False, header = True)
                    counter_nonempty_region +=1
                else:
                    region.to_csv(out_rg, mode='a', index = False, header = False)
            if client:
                print("Shutting down Dask cluster...")
                client.close()

    elif locusbreaker:
        print("Starting LocusBreaker")
        tasks = []
        #Defining the Dask functions for delayed
        traits = pd.read_csv(table)
        @delayed
        def query_spec(uri_path, chrom, trait: str = None, cell: str = None, gene: str = None, type_sumstat:str = "scqtl"):
            with tiledb.open(uri_path, mode="r") as tiledb_data:
                if type_sumstat == "gwas":
                    tiledb_filtered = tiledb_data.query(dims=['CHR','TRAIT','POS'], attrs=['SNPID', 'EAF' , 'BETA', 'SE', 'P', 'N']).df[chrom, trait, :]
                else:
                    tiledb_filtered = tiledb_data.query(dims=['CHR','CELL','GENE','POS'], attrs=['SNPID', 'EAF' , 'BETA', 'SE', 'P', 'N', 'DIST']).df[chrom, cell ,gene , :]
                return tiledb_filtered
            
        @delayed
        def delayed_locus_breaker(tiledb_data, pvalue_sig, pvalue_limit, hole_size, phenovar, category, type_sumstat = "scqtl"):
            # Call locus_breaker with the computed tiledb_data
            return locus_breaker(tiledb_data, pvalue_sig=pvalue_sig, pvalue_limit=pvalue_limit, hole_size=hole_size, phenovar = phenovar, category = category, type_sumstat = type_sumstat)
            
        for ind, row in traits.iterrows():
            chrom = row["CHR"]
            if type_sumstat == "gwas":
                trait = row["TRAIT"]
            else:
                cell,gene = row["TRAIT"].split(":")
                
            task = delayed_locus_breaker(query_spec(uri_path, chrom,trait = trait, cell = cell, gene = gene, type_sumstat = type_sumstat),pvalue_sig=pvalue_sig,pvalue_limit=pvalue_limit,hole_size=hole, phenovar = phenovar, category = category, type_sumstat = type_sumstat)
            tasks.append(task)

        #The batch size to use which is set to the number of workers if Dask is run
        if ctx.obj["workers"]:
            batch_size = ctx.obj["workers"]
        else:
            batch_size = 1
        for i in range(0, len(tasks), batch_size):
                print(f"Batch {i} of {len(tasks)}")
                batch = tasks[i:i+batch_size]
                batch_results = compute(*batch)  # Compute the batch
                for result in batch_results:
                    print(result)
                    #if not len(result) == 0 and not result[0].empty:
                    if result and isinstance(result[0], pd.DataFrame) and not result[0].shape[0] == 0:   
                        interval = result[0]
                        segments = result[1]
                        write_header_interval = not os.path.exists(out_lb + "_interval.csv")
                        write_header_segment = not os.path.exists(out_lb + "_segment.csv")
                        interval.to_csv(out_lb + "_interval.csv", mode="a", index=False, header = write_header_interval)
                        segments.to_csv(out_lb + "_segment.csv", mode="a", index=False, header = write_header_segment)
        if client:
            print("Shutting down Dask cluster...")
            client.close()

    #If no SNP or locusbreker is run only a filtering is done
    else:
        with tiledb.open(uri_path, mode="r") as A:
            if type_sumstat == "gwas":
                tiledb_iterator = A.query(
                    return_incomplete=True,
                    attrs=attr.split(",")
                ).df[chrom, trait , unique_positions]
            else:
                tiledb_iterator = A.query(
                    return_incomplete=True,
                    attrs=attr.split(",")
                ).df[chrom, cell, gene , unique_positions]

            for chunk in tiledb_iterator:
                chunk.to_csv(out_rg + ".csv", mode="a", index=False, header = True)
        print(f"Saved filtered summary statistics in {out_rg}")
        if client:
            print("Shutting down Dask cluster...")
            client.close()
        exit()
    


