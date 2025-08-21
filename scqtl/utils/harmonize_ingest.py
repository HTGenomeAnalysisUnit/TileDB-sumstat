import logging
import json
from pathlib import Path
import os
import pandas as pd
import polars as pl
import numpy as np
import tiledb
import gwaslab as gl

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


class HarmonizationError(Exception):
    """Custom exception for harmonization errors."""
    pass


class Harmonize:
    def __init__(self, mapping_file: str, chunk_size: int, uri: str, type_sumstat: str):
        self.mapping_file = mapping_file
        self.chunk_size = chunk_size
        self.uri = uri
        self.type_sumstat = type_sumstat
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
                tiledb.Attr(name="P", dtype=np.float64, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="N", dtype=np.int64, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
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
                allows_duplicates=True
        )
        tiledb.Array.create(self.uri, schema)
        

    def harmonize(self, file_path, sumstat, trait: str = None, cell: str = None, gene: str = None, N: int = None, qc:bool = False):
        """Load and rename columns, and ensure CHR/POS exist."""
        self.chunk_pl = sumstat.rename(self.mapping_types)

        if not "N" in self.chunk_pl.columns:
            if N is not None:
                self.chunk_pl = self.chunk_pl.with_columns(
                    pl.lit(N).alias("N")
                )
            else:
                raise HarmonizationError("N column is missing and N parameter is not provided")

        # If CHR/POS missing, extract from SNP ID
        if "CHR" not in self.chunk_pl.columns or "POS" not in self.chunk_pl.columns:
            self.chunk_pl = self.chunk_pl.with_columns(
                pl.col("SNPID")
                .str.split_exact(":", 4)
                .struct.rename_fields(["CHR", "POS", "NEA", "EA"])
                .alias("fields")
            ).unnest("fields")
        if "SNPID" not in self.chunk_pl.columns:
            self.chunk_pl = self.chunk_pl.with_columns(
            pl.when(pl.col("NEA") > pl.col("EA"))
            .then(pl.col("BETA"))
            .otherwise(-pl.col("BETA")),
            pl.when(pl.col("NEA") > pl.col("EA"))
            .then(pl.col("EAF"))
            .otherwise(1.0 - pl.col("EAF"))
            )
            self.chunk_pl = self.chunk_pl.with_columns(
                pl.concat_str(
                    pl.lit("chr"),
                    pl.col("CHR")).alias("chr_SNP"))
            self.chunk_pl = self.chunk_pl.with_columns(
                    pl.concat_str(
                    pl.col("chr_SNP"),
                    pl.col("POS"),
                    pl.col("EA"),
                    pl.col("NEA"),
                    separator=":"
                    ).alias("SNPID")
                )

        
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
                "N": np.int64,
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
                "N": np.int64,
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
    
    def qc_sumstat(self, file_path:str):
        directory = self.uri + "_logs"
        filename = os.path.basename(file_path)   # "test.csv.gz"
        # Remove all extensions
        file_name = filename.split('.')[0]       # "test"
    
        if not os.path.isdir(directory):
            os.mkdir(directory)
        sumstat_preqc = self.chunk_pl.to_pandas()
        if self.type_sumstat == "gwas":
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
                 other = ["CELL","GENE","RSID","DIST"])
        sumstat_gl.fix_id()
        sumstat_gl.fix_chr(remove=True)
        sumstat_gl.fix_pos(remove=True)
        sumstat_gl.fix_allele(remove=True)
        sumstat_gl.check_sanity()
        sumstat_gl.check_data_consistency()
        sumstat_gl.remove_dup(mode="m")




        #sumstat_gl.basic_check(n_cores = 4, remove=True, remove_dup=True)

        sumstat_gl.log.save(directory + "/" + file_name)
        self.chunk_pl = pl.from_pandas(sumstat_gl.data)


    def align_alleles(self, pvar_file: str):
        """Align alleles based on pvar file."""
        if not self.pvar_file:
            raise HarmonizationError("pvar_file must be provided to verify alleles order")
        if not Path(self.pvar_file).is_file():
            raise FileNotFoundError(f"pvar_file {self.pvar_file} does not exist")

        pvar_df = pl.read_csv(
            self.pvar_file,
            separator="\t",
            has_header=True,
            new_columns=["CHROM", "POS", "SNPID", "NEA", "EA", "INFO"],
            dtypes={"CHROM": pl.Utf8, "POS": pl.Utf8, "SNPID": pl.Utf8,
                    "REF": pl.Utf8, "ALT": pl.Utf8, "INFO": pl.Utf8},
        )

        self.chunk_pl = self.chunk_pl.join(pvar_df, on="SNPID", how="inner", suffix="_pvar")

        self.chunk_pl = self.chunk_pl.with_columns(
            pl.when(pl.col("NEA") > pl.col("EA"))
            .then(pl.col("BETA"))
            .otherwise(-pl.col("BETA")),
            pl.when(pl.col("NEA") > pl.col("EA"))
            .then(pl.col("EAF"))
            .otherwise(1.0 - pl.col("EAF"))
        )

    def ingest_data(self, file_path):
        """Append harmonized data to TileDB."""
        pl.Config.set_tbl_cols(-1)
        self.chunk_pl = self.chunk_pl.select(self.tiledb_types.keys())
        try:
            tiledb.from_pandas(
                uri=self.uri,
                dataframe=self.chunk_pl.to_pandas(),
                index_dims=self.dimension_tiledb,
                column_types=self.tiledb_types,
                allows_duplicates = False,
                mode="append",
            )
            logger.info(f"Successfully appended chunk to TileDB for file {file_path}")
        except Exception as e:
            logger.error(f"Failed to append chunk to TileDB for file {file_path}: {e}")
            raise

    def create_metadata(self, file_path: str, pvar_file: str = None):
        """Create and store metadata in TileDB."""
        tiledb_existing = tiledb.open(self.uri)
        if "metadata" in tiledb_existing.meta:
            metadata = json.loads(tiledb_existing.meta["metadata"])
        else:
            metadata = {
            "file_path": file_path,
            "type_sumstat": self.type_sumstat,
            "celltype": [],
            "trait": []
            }

        if self.type_sumstat == "qtl":
        # Get unique cell types
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
            
                chr_gene_map = df_cell.group_by("CHR").agg(pl.col("GENE").unique().alias("genes"))
                # Append to metadata, making sure we extend if already exists
                for row in chr_gene_map.iter_rows(named=True):
                    chrom = row["CHR"]
                    genes = row["genes"]
                    if chrom in metadata[cell]:
                        existing = set(metadata[cell][chrom])
                        metadata[cell][chrom].extend([g for g in genes if g not in existing])
                    else:
                        metadata[cell][chrom] = genes   
        else:
            if "TRAIT" not in self.chunk_pl.columns:
                raise HarmonizationError("TRAIT column is missing in the data")
            if self.chunk_pl["TRAIT"].is_empty():
                raise HarmonizationError("TRAIT column is empty in the data")
            traits = self.chunk_pl["TRAIT"].unique().to_list()
            for record in traits:
                if record not in metadata["trait"]:
                    metadata["trait"].append(record)

        with tiledb.open(self.uri, mode='w') as array:
            array.meta["metadata"] = json.dumps(metadata)

        logger.info("Metadata created successfully")