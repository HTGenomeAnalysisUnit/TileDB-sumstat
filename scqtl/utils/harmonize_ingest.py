import logging
import json
from pathlib import Path
import os
import pandas as pd
import polars as pl
import numpy as np
from scipy import stats
import tiledb
import gwaslab as gl
from collections import defaultdict
from scqtl.utils import acat_optimized, compute_pheno_variance


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


class HarmonizationError(Exception):
    """Custom exception for harmonization errors."""
    pass


class Harmonize:
    def __init__(self, mapping_file: str, chunk_size: int, uri: str, type_sumstat: str, pvar_file: str, type_trait: str):
        self.mapping_file = mapping_file
        self.chunk_size = chunk_size
        self.uri = uri
        self.pvar_file = pvar_file
        self.type_sumstat = type_sumstat
        self.type_trait = type_trait
        self.mapping_types = {}
        self.tiledb_types = {}
        self.dimension_tiledb = []
    
    def create_mapping(self):
        df = pd.read_csv(self.mapping_file, header=None, names=["key", "value"])
        if df.empty:
            raise HarmonizationError("Mapping file is empty or not formatted correctly.")
        self.mapping_types = dict(zip(df["key"], df["value"]))
        #check that "BETA", "SE", "AF" are in the vlaues of the mapping_types
        if not all(col in self.mapping_types.values() for col in ["BETA", "SE", "EAF"]):
            raise HarmonizationError("Mapping file must contain BETA, SE, and EAF columns.")
        # Check if CHR and POS or SNP are present
        if not all(col for col in ["CHR", "POS"] if col in self.mapping_types.values()):
            if "SNPID" not in self.mapping_types.values():
                raise HarmonizationError("Mapping file must contain either CHR, and POS or SNPID columns.")
        
    def create_tiledb(self):
        """Create the mapping and dtype definitions."""
        pos_domain = (1, 300000000)  # Example range for genomic positions
        chr_domain = (1, 24)  # Example range for genomic positions
        attrs=[
                tiledb.Attr(name="SNPID", dtype="ascii", filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="RSID", dtype="ascii", filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="EAF", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="BETA", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="SE", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="P", dtype=np.float64, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
            ]
    
        if self.type_sumstat == "gwas":
            
            self.dimension_tiledb = ["CHR", "TRAIT", "POS"]
            dom = tiledb.Domain(
            tiledb.Dim(name="CHR", domain = chr_domain, dtype=np.uint16,  filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            tiledb.Dim(name="TRAIT", dtype="ascii",  var=True, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]) ),
            tiledb.Dim(name="POS", domain = pos_domain, dtype=np.uint32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
            )
        else:
            
            self.dimension_tiledb = ["CHR", "CELL", "GENE", "POS"]
            dom = tiledb.Domain(
                tiledb.Dim(name="CHR", domain = chr_domain, dtype=np.uint16,  filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Dim(name="CELL", dtype="ascii",  var=True, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]) ),
                tiledb.Dim(name="GENE",dtype="ascii", var=True, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]) ),
                tiledb.Dim(name="POS", domain = pos_domain, dtype=np.uint32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
                )
            
            attrs = attrs + [
            tiledb.Attr(name="DIST", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
            ]
        schema = tiledb.ArraySchema(
                domain=dom,
                attrs=attrs,
                sparse=True,
                allows_duplicates=False
        )
        tiledb.Array.create(self.uri, schema)
    
    def align_alleles(self):
        """Align alleles based on pvar file."""
        if not self.pvar_file:
            raise HarmonizationError("pvar_file must be provided to verify alleles order")
        if not Path(self.pvar_file).is_file():
            raise FileNotFoundError(f"pvar_file {self.pvar_file} does not exist")

        pvar_df = pl.read_csv(
            self.pvar_file,
            separator="\t",
            has_header=True,
            dtypes={"CHROM": pl.Utf8, "POS": pl.Utf8, "SNPID": pl.Utf8,
                    "REF": pl.Utf8, "ALT": pl.Utf8},
        )

        self.chunk_pl = self.chunk_pl.join(pvar_df, on="SNPID", how="inner", suffix="_pvar")

        swap = pl.col("ALT")< pl.col("REF")

        self.chunk_pl = self.chunk_pl.with_columns([
            pl.when(swap).then(pl.col("BETA")).otherwise(-pl.col("BETA")),
            pl.when(swap).then(pl.col("EAF")).otherwise(1.0-pl.col("EAF"))
            ])

    def harmonize(self, sumstat, trait: str = None, cell: str = None, gene: str = None, n: int = None, n_controls: int = None,n_cases: int = None, qc:bool = False):
        """Load and rename columns, and ensure CHR/POS exist."""
        self.chunk_pl = sumstat.rename(self.mapping_types)
        
        if self.type_trait == "quant": 
            if not "N" in self.chunk_pl.columns:
                if n is not None:
                    self.chunk_pl = self.chunk_pl.with_columns(
                    pl.lit(n).alias("N")
                    )
                else:
                    raise HarmonizationError("N column is missing and N parameter is not provided")
        elif self.type_trait == "binary": 
            if not "N_CASES" in self.chunk_pl.columns and "N_CONTROLS":
                if n_cases is not None and n_controls is not None:
                    self.chunk_pl = self.chunk_pl.with_columns(
                    pl.lit(float(n_cases)).alias("N_CASES"),
                    pl.lit(float(n_controls)).alias("N_CONTROLS"),
                    pl.lit(float(n_cases) + float(n_controls)).alias("N"),
                    )
                else:
                    raise HarmonizationError("n_cases and n_controls columns are missing and were not provided")
        else:
            raise HarmonizationError("Type of trait must be either binary or quant")




        # If CHR/POS missing, extract from SNP ID
        if "CHR" not in self.chunk_pl.columns or "POS" not in self.chunk_pl.columns:
            self.chunk_pl = self.chunk_pl.with_columns(
                pl.col("SNPID")
                .str.split_exact(":", 4)
                .struct.rename_fields(["CHR", "POS", "A1", "A2"])
                .alias("fields")
            ).unnest("fields")
    
        if "SNPID" in self.chunk_pl.columns:
            self.chunk_pl = self.chunk_pl.drop("SNPID")
                
        swap = pl.col("A1") < pl.col("A2")

        #Start by creating new SNPID aligned
        self.chunk_pl = self.chunk_pl.with_columns([
        pl.when(swap).then(pl.col("A1")).otherwise(pl.col("A2")).alias("EA"),
        # set NEA to the larger allele
        pl.when(swap).then(pl.col("A2")).otherwise(pl.col("A1")).alias("NEA")]
        )

        self.chunk_pl = self.chunk_pl.with_columns(
                pl.concat_str(
                    [
                    pl.col("CHR"),
                    pl.col("POS").cast(pl.Utf8),            # cast POS if numeric
                    pl.col("EA"),       # lexicographically smaller
                    pl.col("NEA")      # lexicographically larger
                    ],
                    separator=":"
                    ).alias("SNPID")
            )

        if self.pvar_file:
            self.align_alleles()
            self.chunk_pl.drop(["REF","ALT"])
        else:
            # flip the sign of BETA when swapping
            self.chunk_pl = self.chunk_pl.with_columns([
                pl.when(swap).then(-pl.col("BETA")).otherwise(pl.col("BETA")).alias("BETA"),
                # flip EAF to 1 - EAF when swapping
                pl.when(swap).then(1.0 - pl.col("EAF")).otherwise(pl.col("EAF")).alias("EAF"),
                # set EA to the smaller allele
            ])
        self.chunk_pl = self.chunk_pl.drop(["A1","A2"])

        self.chunk_pl = self.chunk_pl.with_columns(
                    pl.concat_str(
                        pl.lit("chr"),
                        pl.col("SNPID")).alias("SNPID"))
    
        if self.type_sumstat=="gwas":
            self.tiledb_types = {
                    "CHR": np.uint16,
                    "TRAIT": str,
                    "POS": np.uint32,
                    "SNPID": str,
                    "RSID": str,
                    "EAF": np.float32,
                    "BETA": np.float32,
                    "SE": np.float32,
                    "P": np.float64,
                }
            if "TRAIT" not in self.chunk_pl.columns:
                self.chunk_pl = self.chunk_pl.with_columns(
                    pl.lit(trait).alias("TRAIT")
                )
        else:
            self.tiledb_types = {
                "CHR": np.uint16,
                "CELL": str,
                "GENE": str,
                "POS": np.uint32,
                "SNPID": str,
                "RSID": str,
                "DIST": np.int64,
                "EAF": np.float32,
                "BETA": np.float32,
                "SE": np.float32,
                "P": np.float64,
            }
            if "CELL" not in self.chunk_pl.columns:
                self.chunk_pl = self.chunk_pl.with_columns(
                    pl.lit(cell).alias("CELL")
                )
            if "GENE" not in self.chunk_pl.columns:
                self.chunk_pl = self.chunk_pl.with_columns(
                    pl.lit(gene).alias("GENE")
                )
        if "RSID" not in self.chunk_pl.columns:
                self.chunk_pl = self.chunk_pl.with_columns(
                    pl.lit("None").alias("RSID")
                )
        if "LOG10P" in self.chunk_pl.columns:
            self.chunk_pl = self.chunk_pl.with_columns(
                    (10 ** (-pl.col("LOG10P"))).alias("P")
                )
        
        #Calculate p-value from z-score
        self.chunk_pl = self.chunk_pl.drop('P')
        self.chunk_pl = self.chunk_pl.with_columns(
            (pl.col("BETA") / pl.col("SE")).pow(2).map_batches(
            lambda x: pl.Series(stats.chi2.sf(x.to_numpy(), df=1)),
            return_dtype=pl.Float64
            ).alias('P')
            )
    
    def qc_sumstat(self, file_path:str):
        directory = self.uri + "_logs"
        filename = os.path.basename(file_path)   # "test.csv.gz"
        # Remove all extensions
        file_name = filename.split('.')[0]       # "test"
    
        if not os.path.isdir(directory):
            os.mkdir(directory)
        sumstat_preqc = self.chunk_pl.to_pandas()
        if self.type_sumstat == "gwas":
            if self.type_trait== "quant":
                sumstat_gl =gl.Sumstats(sumstat_preqc,
                    snpid="SNPID",
                    chrom="CHR",
                    pos="POS",
                    eaf="EAF",
                    beta="BETA",
                    se="SE",
                    p="P",
                    n="N",
                    ea = "EA",
                    nea = "NEA",
                    other = ["TRAIT","RSID"])
            else:
                sumstat_gl =gl.Sumstats(sumstat_preqc,
                    snpid="SNPID",
                    chrom="CHR",
                    pos="POS",
                    eaf="EAF",
                    beta="BETA",
                    se="SE",
                    p="P",
                    n="N",
                    ncase = "N_CASES",
                    ncontrol = "N_CONTROLS",
                    ea = "EA",
                    nea = "NEA",
                    other = ["TRAIT","RSID"])

        else:
            sumstat_gl =gl.Sumstats(sumstat_preqc,
                 snpid="SNPID",
                 chrom="CHR",
                 pos="POS",
                 eaf="EAF",
                 beta="BETA",
                 se="SE",
                 p="P",
                 n="N",
                 other = ["CELL","GENE","RSID","DIST"])
        #sumstat_gl.fix_id()
        sumstat_gl.fix_chr(remove=True)
        sumstat_gl.fix_pos(remove=True)
        sumstat_gl.fix_allele(remove=True)
        sumstat_gl.check_sanity()
        sumstat_gl.check_data_consistency()
        #sumstat_gl.remove_dup(mode="m")
        #sumstat_gl.basic_check(n_cores = 4, remove=True, remove_dup=True)

        sumstat_gl.log.save(directory + "/" + file_name)
        self.chunk_pl = pl.from_pandas(sumstat_gl.data)
    
    def ingest_data(self, file_path):
        """Append harmonized data to TileDB."""
        pl.Config.set_tbl_cols(-1)
        if self.type_sumstat == "gwas":
            dedup_keys = ["CHR", "POS", "TRAIT"]
        else:
            dedup_keys = ["CHR", "POS", "CELL", "GENE"]

        # Option A (recommended): window count -> keep only rows whose group count == 1
        self.chunk_pl = (
            self.chunk_pl
                .with_columns(pl.count().over(dedup_keys).alias("_grp_count"))
                .filter(pl.col("_grp_count") == 1)
                .drop("_grp_count")
                )
        chunk_pl_ingest = self.chunk_pl.select(self.tiledb_types.keys())
        try:
            tiledb.from_pandas(
                uri=self.uri,
                dataframe=chunk_pl_ingest.to_pandas(),
                index_dims=self.dimension_tiledb,
                column_types=self.tiledb_types,
                allows_duplicates = False,
                mode="append"
            )
            logger.info(f"Successfully appended chunk to TileDB for file {file_path}")
        except Exception as e:
            logger.error(f"Failed to append chunk to TileDB for file {file_path}: {e}")
            raise

    def create_metadata(self, file_path: str):
        """Create and store metadata in TileDB."""
        #tiledb_existing = tiledb.open(self.uri)
        #if "metadata" in tiledb_existing.meta:
        #    metadata = json.loads(tiledb_existing.meta["metadata"])
        #else:
        metadata = {
            "file_path": file_path,
            "CELL": [],
            "trait": []
            }
        if self.type_sumstat == "qtl":
        # Get unique cell types

            self.chunk_pl = self.chunk_pl.group_by(["CHR" ,"CELL", "GENE"]).agg(
                pl.col("P").map_batches(
                    lambda s: pl.Series([acat_optimized(s)]),
                    return_dtype=pl.Float64
                ).alias("ACAT_P"),
                pl.col("N").first().alias("N")
                )
            self.chunk_pl = self.chunk_pl.with_columns(
                self.chunk_pl["ACAT_P"].list.first().alias("ACAT_P_scalar")
            )

            celltypes = self.chunk_pl["CELL"].unique().to_list()
            if not celltypes:
                raise HarmonizationError("No cell types found in the data")

            # Loop over each cell type
            for cell in celltypes:
                if cell not in metadata["celltype"]:
                    metadata["celltype"].append(cell)
                    metadata[cell] = {}

                # Filter by this cell type
                df_cell = self.chunk_pl.filter(self.chunk_pl["CELL"] == cell)
                # Group by chromosome and collect unique genes
                if self.type_trait == "quant":
                    chr_gene_map = df_cell.group_by(["CHR", "GENE"]).agg([
                        pl.col("ACAT_P_scalar").first().alias("ACAT"),
                        pl.col("N").first()
                        ])
                    for row in chr_gene_map.iter_rows(named=True):
                        chrom = row["CHR"]
                        gene = row["GENE"]
                        acat = row["ACAT"]
                        n = float(row["N"])
                        if chrom not in metadata[cell]:
                            metadata[cell][chrom] = {}
                            metadata[cell][chrom][gene] = {
                                "N": n,
                                "ACAT": acat,
                                "PHENO_VAR": 1
                                }
                else:
                    chr_gene_map = df_cell.group_by(["CHR", "GENE"]).agg([
                        pl.col("ACAT_P_scalar").first().alias("ACAT"),
                        pl.col("N").first(),
                        pl.col("N_CASES").first(),
                        pl.col("N_CONTROLS").first()
                        ])
                    for row in chr_gene_map.iter_rows(named=True):
                        chrom = row["CHR"]
                        gene = row["GENE"]
                        acat = row["ACAT"]
                        n = float(row["N"])
                        ncases = float(row["N_CASES"])
                        ncontrols = float(row["N_CONTROLS"])
                        if chrom not in metadata[cell]:
                            metadata[cell][chrom] = {}
                            metadata[cell][chrom][gene] = {
                                "N": n,
                                "N_CASES":ncases,
                                "N_CONTROLS":ncontrols,
                                "ACAT": acat,
                                "PHENO_VAR": 1
                                }                    
        else:
            if "TRAIT" not in self.chunk_pl.columns:
                raise HarmonizationError("TRAIT column is missing in the data")
            if self.chunk_pl["TRAIT"].is_empty():
                raise HarmonizationError("TRAIT column is empty in the data")
            traits = self.chunk_pl["TRAIT"].unique().to_list()
            for trait in traits:
                if trait not in metadata["trait"]:
                    metadata["trait"].append(trait)
                    df_trait = self.chunk_pl.filter(self.chunk_pl["TRAIT"] == trait)
                    pheno_var = compute_pheno_variance(df_trait, self.type_trait)
                    n = df_trait["N"].unique().to_list()[0]
                    metadata[trait] = {"N":n, "PHENO_VAR":pheno_var}
                    if self.type_trait == "binary":
                        n_cases = df_trait["N_CASES"].unique().to_list()[0]
                        n_controls = df_trait["N_CONTROLS"].unique().to_list()[0]
                        metadata[trait]["N_CASES"] = n_cases
                        metadata[trait]["N_CONTROLS"] = n_controls
                    
        f = open(f'{self.uri}_metadata.json', 'a')

        with open(f'{self.uri}_metadata.json', 'a') as f:
            json.dump(metadata, f)

        logger.info("Metadata created successfully")
    
    def ingest_metadata(self):
        metadata_records = []
        with open(f'{self.uri}_metadata.json') as f:
            text = f.read()
        
        tdb = tiledb.open(self.uri, 'w')
        decoder = json.JSONDecoder()
        idx = 0
        while idx < len(text):
            obj, idx = decoder.raw_decode(text, idx)
            metadata_records.append(obj)

        # Create a merged dictionary
        merged_metadata = {
            "file_path": [],
            "trait": [],
            "celltype": []
        }

        # For storing per-cell/chromosome info if present
        cell_chrom_data = defaultdict(lambda: defaultdict(list))
        
        for record in metadata_records:
            # Merge file paths
            if "file_path" in record and record["file_path"] not in merged_metadata["file_path"]:
                merged_metadata["file_path"].append(record["file_path"])
            # Merge traits
            if "trait" in record and len(record["trait"]) > 0:
                for t in record["trait"]:
                    if t not in merged_metadata["trait"]:
                        merged_metadata["trait"].append(t)
                    if t in record:
                        merged_metadata[t] = record[t]
            # Merge celltypes and their per-chromosome data
            if "cell" in record and len(record["celltype"]) > 0:
                for cell in record["celltype"]:
                    if cell not in merged_metadata["celltype"]:
                        merged_metadata["celltype"].append(cell)
                    # Merge chromosome-level info if exists
                    if cell in record:
                        for chrom, genes in record[cell].items():
                            cell_chrom_data[cell][chrom] = {}
                            for g in genes:
                                cell_chrom_data[cell][chrom][g] = genes[g]

        # Add per-cell/chromosome info
        for cell, chrom_dict in cell_chrom_data.items():
            merged_metadata[cell] = chrom_dict

        rows = []
        tdb = tiledb.open(self.uri, 'w')
        tdb.meta["metadata"] = json.dumps(merged_metadata)
        # Loop over all celltypes in metadata
        for celltype in merged_metadata["celltype"]:
            if celltype in merged_metadata:  # make sure the key exists in dict
                groups = merged_metadata[celltype]
                for group, values in groups.items():  # e.g. group "11"
                    for gene_id in values:
                        if self.type_trait=="quant":
                            rows.append({
                                "CHR": group,
                                "CELL": celltype,
                                "GENE": gene_id,
                                "N": values[gene_id]["N"],
                                "ACAT": values[gene_id]["ACAT"],
                                "PHENOVAR": values[gene_id]["PHENO_VAR"],
                            })
                        else:
                            rows.append({
                                "CHR": group,
                                "CELL": celltype,
                                "GENE": gene_id,
                                "N": values[gene_id]["N"],
                                "N_CASES": values[gene_id]["N_CASES"],
                                "N_CONTROLS": values[gene_id]["N_CONTROLS"],
                                "ACAT": values[gene_id]["ACAT"],
                                "PHENOVAR": values[gene_id]["PHENO_VAR"],
                            })


        df = pd.DataFrame(rows)
        df.to_csv(f"{self.uri}_metadata.csv",index = False)
        

