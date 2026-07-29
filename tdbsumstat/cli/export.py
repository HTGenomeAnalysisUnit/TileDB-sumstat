import tiledb
import click
import cloup
import pandas as pd
from tdbsumstat.utils.locusbreaker_plpl import locusbreaker_plpl
from tdbsumstat.utils import acat_optimized
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
    cloup.option("--uri-path", default = None, type=str, help = "path of TileDB"),
    cloup.option("--table-regions", default = None, type=str, help = "Regions to interrogate from a table"),
    cloup.option("--trait-list", default = None, type=str, help = "List of entire traits to filter"),
    cloup.option("--attr", default = "P,SNPID,EAF,BETA,SE", type=str, help = "Attributes to output"),
    cloup.option("--export-meta", is_flag = True, default = False, type=str, help = "Get metadata from TileDB"),
    cloup.option("--mac", default = 0, type=int, help = "Filter for MAC when recomputing metadata"),
    cloup.option("--recompute-meta", is_flag = True, default = False, type=str, help = "Recompute metadata after applying filters (Does not modify data within the TileDB)"),
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
        uri_path: str,
        type_sumstat: str,
        table_regions: str,
        trait_list:str,
        mac:int,
        attr: str,
        snp: str,
        export_meta:bool,
        recompute_meta:bool,
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
    tiledb_export = tiledb.open(uri_path, mode="r")
    metadata = json.loads(tiledb_export.meta["merged_metadata"])
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
        for trait in metadata["traits"]:
            rows.append({
                    "TRAIT": trait,
                    **metadata[trait]
                    })
        df_meta = pl.DataFrame(rows)
    if snp: 
        snp_list = pd.read_csv(snp, dtype = {"CHR":int, "POS":np.uint32, "TRAIT":str})        
        if type_sumstat == "gwas":
            trait_list = snp_list['TRAIT'].unique().tolist()
        else:
            snp_list[["CELL","GENE"]] = snp_list["TRAIT"].str.split(":", n = 2, expand = True)
            trait_list = snp_list['CELL'].unique().tolist()
        for trait in trait_list:
            if type_sumstat == "gwas":
                chrom_list = snp_list[snp_list['TRAIT']==trait]['CHR'].unique().tolist()
            else:
                chrom_list = snp_list[snp_list['CELL']==trait]['CHR'].unique().tolist()
            for chrom in chrom_list:
                if type_sumstat == "gwas":
                    snp_list_refined = snp_list[(snp_list['CHR']==chrom) & (snp_list['TRAIT']==trait)]['POS'].unique().tolist()
                    tiledb_query = tiledb_export.query(
                        	attrs=attr.split(","),
                            return_arrow = True
                    	).df[chrom, trait , :]
                    tiledb_query_pl = pl.from_arrow(tiledb_query)
                    tiledb_query_pd = tiledb_query_pl.filter(pl.col("POS").is_in(snp_list_refined)).to_pandas()
                else:
                    snp_list_refined = snp_list[(snp_list['CHR']==chrom) & (snp_list['CELL']==trait)]['POS'].unique().tolist()
                    gene_list = snp_list[(snp_list['CHR']==chrom) & (snp_list['CELL']==trait)]['GENE'].unique().tolist()
                    tiledb_query = tiledb_export.query(
                        attrs=attr.split(",")
                    ).df[chrom, trait , gene_list, :]
                    tiledb_query_pl = pl.from_arrow(tiledb_query)
                    tiledb_query_pd = tiledb_query_pl.filter(pl.col("POS").is_in(snp_list_refined)).to_pandas()
                tiledb_query_pd.to_csv(f"{out}_{trait}_{chrom}.csv", mode = "a", index=False)
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
        def query_spec(uri_path, chrom:int, trait: str = None, cell: str = None, gene: str = None, type_sumstat:str = "scqtl", region = slice(None)):
            with tiledb.open(uri_path, mode="r") as tiledb_data:
                if type_sumstat == "gwas":
                    tiledb_filtered = tiledb_data.query(dims=['CHR','TRAIT','POS']).df[chrom, trait, region]
                else:
                    tiledb_filtered = tiledb_data.query(dims=['CHR','CELL','GENE','POS'], return_arrow=True).df[chrom, cell ,gene , region]
                return tiledb_filtered
            #return locusbreaker_plpl(tiledb_data, maf = maf, pvalue_sig=pvalue_sig, pvalue_limit=pvalue_limit, locus_max_size = locus_max_size, 
            #                         hole_size=hole_size, cis_trans_lb = cis_trans_lb, type_sumstat = type_sumstat, metadata = metadata)
        traits = pd.read_csv(table_lb)
        traits = traits.astype({'CHR': 'int16'})
        for index, trait in traits.iterrows():
                if "SIG" in traits.columns:
                    pvalue_sig = trait["SIG"]
                    pvalue_limit = trait["LIM"]
                region = slice(None)
                if "START" in traits.columns and "END" in traits.columns:
                    region=slice(trait["START"],trait["END"])
                if type_sumstat == "gwas":
                    query = query_spec(uri_path, trait["CHR"], trait=trait["TRAIT"], type_sumstat=type_sumstat, region=region)
                else:
                    cell,genes = trait["TRAIT"].split(":")
                    query = query_spec(uri_path, trait["CHR"], cell=cell, gene=genes, type_sumstat=type_sumstat, region=region)

                result = locusbreaker_plpl(query,
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
    elif export_meta:
        tiledb_db = tiledb.open(uri_path, mode="r")
        tiledb_meta = json.loads(tiledb_db.meta['merged_metadata'])
        if type_sumstat == "qtl":
            rows = []
            for cell_type in tiledb_meta["CELL"]:
                samples = tiledb_meta.get(cell_type, {})
                for sample_id, genes in samples.items():
                    for gene_id, metrics in genes.items():
                        rows.append({
                            "CELL": cell_type,
                            "CHR": sample_id,
                            "GENE": gene_id,
                        **metrics
                        })
        else:
            for trait in tiledb_meta["traits"]:
                samples = tiledb_meta.get(trait, {})
                for chrom, genes in samples.items():
                    for gene_id, metrics in genes.items():
                        rows.append({
                            "CELL": cell_type,
                            "CRH": chrom,
                        **metrics
                        })   

        df = pd.DataFrame(rows)
        df.to_csv(f"{out}_meta.csv", index = False)
    
    elif recompute_meta:
        trait_list_pd = pd.read_csv(trait_list)
        for record, trait in trait_list_pd.iterrows():
            if type_sumstat == "gwas":
                tiledb_query = tiledb_export.query().df[int(trait['CHR']), trait['TRAIT'].to_string() , :]
                n = df_meta.filter(pl.col('TRAIT')==trait['TRAIT']).select('N')["N"][0]
            else:
                #cell,gene = trait['TRAIT'].split(':')
                print(trait['CELL'])
                tiledb_query = tiledb_export.query().df[int(trait['CHR']), trait['CELL'], : , :]
                
                n = df_meta.filter(pl.col('CELL')==trait['CELL']).select('N')["N"][0]
                print(n)
            if "N" not in tiledb_query.columns:
                tiledb_query["N"] = n
                print(tiledb_query)
            tiledb_query_pl = pl.from_pandas(tiledb_query)
            
            tiledb_query_pl= tiledb_query_pl.with_columns(
                (2 * pl.col("N") * pl.min_horizontal("EAF", (1 - pl.col("EAF"))))
                .alias("MAC")
                 ).filter(pl.col("MAC") > mac)
            if type_sumstat == "gwas":
                chr_gene_agg = tiledb_query_pl.group_by(["CHR","TRAIT"]).agg([
                    pl.col("P").map_batches(
                        lambda s: pl.Series([acat_optimized(s)]),
                        return_dtype=pl.Float64
                    ).alias("ACAT_LIST"),
                    pl.col("N").first().alias("N"),
                    pl.min("P").alias("MIN_P")
                ])
                chr_gene_agg = chr_gene_agg.with_columns(
                    pl.col("ACAT_LIST").list.first().alias("ACAT")
                )
            else:
                chr_gene_agg = tiledb_query_pl.group_by(["CHR","CELL","GENE"]).agg([
                    pl.col("P").map_batches(
                        lambda s: pl.Series([acat_optimized(s)]),
                        return_dtype=pl.Float64
                    ).alias("ACAT_LIST"),
                    pl.col("N").first().alias("N"),
                    pl.min("P").alias("MIN_P")
                ])
                chr_gene_agg = chr_gene_agg.with_columns(
                    pl.col("ACAT_LIST").list.first().alias("ACAT"),
                ).drop("ACAT_LIST")

            chr_gene_agg_pd = chr_gene_agg.to_pandas()
            chr_gene_agg_pd.to_csv(f"{out}_batch_{batch_name}_metadata.csv", mode="a", index=False)
    else:
        trait_list_pd = pd.read_csv(trait_list)
        with tiledb.open(uri_path, mode="r") as A:
            if type_sumstat == "gwas":
                trait_list_np = trait_list_pd["TRAIT"].to_list()
                tiledb_iterator = A.query(
                    return_incomplete=True,
                    attrs=attr.split(",")
                ).df[:, trait_list_np , :]
            else:
                trait_list_pd[['cell','gene']] = trait_list_pd['TRAIT'].str.split('~', expand = True)
                cells = trait_list_pd['cell'].to_list()
                gene = trait_list_pd['gene'].to_list()
                tiledb_iterator = A.query(
                    return_incomplete=True,
                    attrs=attr.split(",")
                ).df[:, cells, gene , :]

            for chunk in tiledb_iterator:
                chunk.to_csv(f"{out}_{batch_name}.csv", mode="a", index=False, header = False)
        print(f"Saved filtered summary statistics in {out}")
