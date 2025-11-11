import tiledb
import click
import cloup
import pandas as pd
from scqtl.utils.locusbreaker_plpl import locusbreaker_plpl
import os
import json
import polars as pl
import random
import numpy as np

help_doc = """
Query TileDB database and export data.
"""

@cloup.command("export", no_args_is_help=True, help=help_doc)
@cloup.option_group(
    "Options for querying specific chromosomes, cells, genes or positions in the TileDB",
    cloup.option("--tiledb-path", default = None, type=str, help = "path of TileDB"),
    cloup.option("--table-regions", default = None, type=str, help = "Regions to interrogate from a table"),
    cloup.option("--trait", default = None, type=str, help = "Trait to filter the for an entire summary statistics"),
    cloup.option("--attr", default = "P,SNPID,EAF,BETA,SE", type=str, help = "Attributes to output"),
    cloup.option("--snp", default = None, type=str, help = "List of SNPs to interrogate taken from a txt file. Please check README for details on the format of this file"),
    cloup.option("--batch-name", default = None, type=str, help = "Name of the batch")
)
@cloup.option_group(
    "Options for Locusbreaker",
    cloup.option("--locusbreaker", is_flag=True, type=bool, default = False, help="Option to run locusbreaker"),
    cloup.option("--hole-lb", default = 250000, type=int, help = "Minimum pair-base distance between SNPs in different loci (default: 250000)"),
    cloup.option("--maf-lb", default = 0.001, type=float, help = "The MAF to filter the TILEDB before locusbreaker is run"),
    cloup.option("--locus-max-size-lb", default = 3000000, type=float, help = "The maximum size allowed for the locus. Default: 1Mb"),
    cloup.option("--cis-trans-lb",default = "cis",type=str,  help = "If locusbreaker run on cis or trans QLTs"),
    cloup.option("--table-lb", default = None, type=str, help = "Path of the table to provide"),
    cloup.option("--type-sumstat", default = None, type=str, help = "Type of summary data")
)
@cloup.option_group(
    "Options for output",
    cloup.option("--out", default = "out", type=str, help = "Output path with file name where results will be stored"),
)

@click.pass_context
def export(
        ctx,
        tiledb_path: str,
        type_sumstat: str,
        table_regions: str,
        trait:str,
        attr: str,
        snp: str,
        locusbreaker: bool,
        maf_lb: float,
        cis_trans_lb: str,
        table_lb: str,
        hole_lb: int,
        out: str,
        locus_max_size_lb: int,
        batch_name:str
        ):
    #Open connection with TileDB
    tiledb_export = tiledb.open(tiledb_path, mode="r")
    metadata = json.loads(tiledb_export.meta["metadata"])
    rows = []
    if type_sumstat == "qtl":
        for cell in metadata["CELL"]:
            for chrom, genes in metadata[cell].items():
                for gene, stats in genes.items():
                    rows.append({
                        "CELL": cell,
                        "CHR": chrom,
                        "GENE": gene,
                        **stats
                    })
        df_meta = pl.DataFrame(rows)  
        df_meta = df_meta.with_columns(pl.col("CHR").cast(pl.UInt16))
    else:
        for trait in metadata["trait"]:
            rows.append({
                    "TRAIT": trait,
                    **metadata[trait]
                    })
        df_meta = pl.DataFrame(rows)
    #Print only the schema of the tiledb
    #Get list of genes, cell type and positions or create ones
    unique_positions = slice(None)
    if not trait:
        trait = slice(None)
    #Intersect the tiledb with a list of SNPs
    if snp: 
        snp_list = pd.read_csv(snp, dtype = {"CHR":int, "POS":np.uint32, "TRAIT":str})        
        if type_sumstat == "gwas":
            header_file = "CHR,TRAIT,POS,SNPID" + "," + attr
            trait_list = snp_list['TRAIT'].unique().tolist()
        else:
            header_file = "CHR,POS,SNPID" + "," + attr + ",TRAIT"
            snp_list[["CELL","GENE"]] = snp_list["TRAIT"].str.split(":", n = 2, expand = True)
            trait_list = snp_list['CELL'].unique().tolist()
        out_file = out + ".csv"
        with open(out_file, "w") as f:
            f.write(header_file + "\n")
        for trait in trait_list:
            if type_sumstat == "gwas":
                chrom_list = snp_list[snp_list['TRAIT']==trait]['CHR'].unique().tolist()
            else:
                chrom_list = snp_list[snp_list['CELL']==trait]['CHR'].unique().tolist()
            for chrom in chrom_list:
                if type_sumstat == "gwas":
                    snp_list_refined = snp_list[(snp_list['CHR']==chrom) & (snp_list['TRAIT']==trait)]['POS'].unique().tolist()
                else:
                    snp_list_refined = snp_list[(snp_list['CHR']==chrom) & (snp_list['TRAIT']==trait)]['POS'].unique().tolist()
                if type_sumstat == "gwas":
                    tiledb_query = tiledb_export.query(
                        	attrs=attr.split(",")
                    	).df[chrom, trait ,snp_list_refined]
                else:
                    tiledb_query = tiledb_export.query(
                        attrs=attr.split(",")
                    ).df[chrom, trait , :, snp_list_refined]
                    tiledb_query["TRAIT"] = tiledb_query["CELL"] + ":" + tiledb_query["GENE"]
                    tiledb_query = tiledb_query.drop(['CELL', 'GENE'], axis = 1)
                merged_df = tiledb_query.merge(snp_list.drop(['CELL','GENE'], axis = 1), on = ["CHR","TRAIT", "POS"])
                if not batch_name:
                        batch_name = random.randint(1, 10000000)
                merged_df.to_csv(f"{out}_batch_{batch_name}.csv", mode="a", index=False, header = False)
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
                    region.to_csv(out, mode='a', index = False, header = True)
                    counter_nonempty_region +=1
                else:
                    region.to_csv(out, mode='a', index = False, header = False)
    elif locusbreaker:
        print("Starting LocusBreaker")
        def query_spec(tiledb_path, chrom:int, trait: str = None, cell: str = None, gene: str = None, type_sumstat:str = "scqtl"):
            with tiledb.open(tiledb_path, mode="r") as tiledb_data:
                if type_sumstat == "gwas":
                    tiledb_filtered = tiledb_data.query(dims=['CHR','TRAIT','POS']).df[chrom, trait, :]
                else:
                    tiledb_filtered = tiledb_data.query(dims=['CHR','CELL','GENE','POS'], return_arrow=True).df[chrom, cell ,gene , :]
                return tiledb_filtered
            
        def locus_breaker(tiledb_data, maf, pvalue_sig, pvalue_limit, locus_max_size, hole_size, cis_trans_lb, metadata, type_sumstat = "scqtl"):
            return locusbreaker_plpl(tiledb_data, maf = maf, pvalue_sig=pvalue_sig, pvalue_limit=pvalue_limit, locus_max_size = locus_max_size, 
                                     hole_size=hole_size, cis_trans_lb = cis_trans_lb, type_sumstat = type_sumstat, metadata = metadata)
        traits = pd.read_csv(table_lb)
        traits = traits.astype({'CHR': 'int16'})
        for index, trait in traits.iterrows():
                if "SIG" in traits.columns:
                    pvalue_sig = trait["SIG"]
                    pvalue_limit = trait["LIM"]

                if type_sumstat == "gwas":
                    query = query_spec(tiledb_path, trait["CHR"], trait=trait["TRAIT"], type_sumstat=type_sumstat)
                else:
                    cell,genes = trait["TRAIT"].split(":")
                    query = query_spec(tiledb_path, trait["CHR"], cell=cell, gene=genes, type_sumstat=type_sumstat)

                result = locus_breaker(query,
                                     maf=maf_lb, 
                                     pvalue_sig=pvalue_sig,
                                     pvalue_limit=pvalue_limit, 
                                     locus_max_size=locus_max_size_lb, 
                                     hole_size=hole_lb, 
                                     cis_trans_lb=cis_trans_lb, 
                                     type_sumstat=type_sumstat,
                                     metadata = df_meta)
        
                if not len(result) == 0 and not result[0].empty:
                    if result and isinstance(result[0], pd.DataFrame) and not result[0].shape[0] == 0:   
                            interval = result[0]
                            segments = result[1]
                            if not batch_name:
                                batch_name = random.randint(1, 10000000)

                            write_header_interval = not os.path.exists(f"{out}_batch_{batch_name}_interval.csv")
                            write_header_segment = not os.path.exists(f"{out}_batch_{batch_name}_segment.csv")
                            interval.to_csv(f"{out}_batch_{batch_name}_interval.csv", mode="a", index=False, header = write_header_interval)
                            segments.to_csv(f"{out}_batch_{batch_name}_segment.csv", mode="a", index=False, header = write_header_segment)
    else:
        with tiledb.open(tiledb_path, mode="r") as A:
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
                chunk.to_csv(out + ".csv", mode="a", index=False, header = True)
        print(f"Saved filtered summary statistics in {out}")
