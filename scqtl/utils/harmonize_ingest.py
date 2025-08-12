import logging
import json
from pathlib import Path

import pandas as pd
import polars as pl
import numpy as np
import tiledb

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
        self.chunk_pl = None  # will hold the harmonized dataframe
    
    def create_mapping(self):
        df = pd.read_csv(self.mapping_file, header=None, names=["key", "value"])
        if df.empty:
            raise HarmonizationError("Mapping file is empty or not formatted correctly.")
        self.mapping_types = dict(zip(df["key"], df["value"]))
        print(f"Mapping types: {self.mapping_types}")
        print(f"Mapping values: {self.mapping_types.values()}")
        #check that "BETA", "SE", "AF" are in the vlaues of the mapping_types
        if not all(col in self.mapping_types.values() for col in ["BETA", "SE", "AF"]):
            raise HarmonizationError("Mapping file must contain BETA, SE, and AF columns.")
        # Check if CHR and POS or SNP are present
        if not all(col for col in ["CHR", "POS"] if col in self.mapping_types.values()):
            if "SNP" not in self.mapping_types.values():
                raise HarmonizationError("Mapping file must contain either CHR, and POS or SNP columns.")
        
    def create_tiledb(self):
        """Create the mapping and dtype definitions."""
        pos_domain = (1, 300000000)  # Example range for genomic positions
        chr_domain = (1, 24)  # Example range for genomic positions
        attrs=[
                tiledb.Attr(name="SNP", dtype="ascii", filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="RSID", dtype="ascii", filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="AF", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="BETA", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="SE", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="P", dtype=np.float64, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="N", dtype=np.int64, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
            ]
    
        if self.type_sumstat == "gwas":
            self.tiledb_types = {
                "CHR": np.uint16,
                "TRAIT": str,
                "POS": np.uint32,
                "SNP": str,
                "RSID": str,
                "AF": np.float32,
                "BETA": np.float32,
                "SE": np.float32,
                "P": np.float64,
                "N": np.int64,
            }
            self.dimension_tiledb = ["CHR", "TRAIT", "POS"]
            dom = tiledb.Domain(
            tiledb.Dim(name="CHR", domain = chr_domain, dtype=np.uint16,  filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            tiledb.Dim(name="TRAIT", dtype="ascii",  var=True, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]) ),
            tiledb.Dim(name="POS", domain = pos_domain, dtype=np.uint32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
            )
        else:
            self.tiledb_types = {
                "CHR": np.uint16,
                "CELL": str,
                "GENE": str,
                "POS": np.uint32,
                "SNP": str,
                "RSID": str,
                "DIST": np.int64,
                "AF": np.float32,
                "BETA": np.float32,
                "SE": np.float32,
                "P": np.float64,
                "N": np.int64,
            }
            self.dimension_tiledb = ["CHR", "CELL", "GENE", "POS"]
            dom = tiledb.Domain(
                tiledb.Dim(name="CHR", domain = chr_domain, dtype=np.uint16,  filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Dim(name="CELL", dtype="ascii",  var=True, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]) ),
                tiledb.Dim(name="GENE",dtype="ascii", var=True, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]) ),
                tiledb.Dim(name="POS", domain = pos_domain, dtype=np.uint32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
                )
            
            attributes = attributes + [
            tiledb.Attr(name="DIST", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
            ]
        schema = tiledb.ArraySchema(
                domain=dom,
                attrs=attrs,
                sparse=True,
                allows_duplicates=True
        )
        tiledb.Array.create(self.uri, schema)
        

    def harmonize(self, file_path, trait: str = None, cell: str = None, gene: str = None, N: int = None):
        """Load and rename columns, and ensure CHR/POS exist."""
        self.chunk_pl = pl.read_csv(
            file_path,
            separator="\t",
            columns=list(self.mapping_types.keys()),
            low_memory=True,
        )
        self.chunk_pl = self.chunk_pl.rename(self.mapping_types)

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
                pl.col("SNP")
                .str.split_exact(":", 4)
                .struct.rename_fields(["CHR", "POS", "REF", "ALT"])
                .alias("fields")
            ).unnest("fields")
        if "SNP" not in self.chunk_pl.columns:
            self.chunk_pl = self.chunk_pl.with_columns(
            pl.when(pl.col("REF") > pl.col("ALT"))
            .then(pl.col("BETA"))
            .otherwise(-pl.col("BETA")),
            pl.when(pl.col("REF") > pl.col("ALT"))
            .then(pl.col("AF"))
            .otherwise(1.0 - pl.col("AF"))
            )
            self.chunk_pl = self.chunk_pl.with_columns(
                pl.concat_str(
                    pl.lit("chr"),
                    pl.col("CHR")).alias("chr_SNP"))
            self.chunk_pl = self.chunk_pl.with_columns(
                    pl.concat_str(
                    pl.col("chr_SNP"),
                    pl.col("POS"),
                    pl.col("ALT"),
                    pl.col("REF"),
                    separator=":"
                    ).alias("SNP")
                )

        
        if self.type_sumstat=="gwas":
            if "TRAIT" not in self.chunk_pl.columns:
                self.chunk_pl = self.chunk_pl.with_columns(
                    pl.lit(trait).alias("TRAIT")
                )
        else:
            if "CELL" not in self.chunk_pl.columns:
                self.chunk_pl = self.chunk_pl.with_columns(
                    pl.lit(cell).alias("CELL")
                )
            if "GENE" not in self.chunk_pl.columns:
                self.chunk_pl = self.chunk_pl.with_columns(
                    pl.lit(gene).alias("GENE")
                )
        

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
            new_columns=["CHROM", "POS", "SNP", "REF", "ALT", "INFO"],
            dtypes={"CHROM": pl.Utf8, "POS": pl.Utf8, "SNP": pl.Utf8,
                    "REF": pl.Utf8, "ALT": pl.Utf8, "INFO": pl.Utf8},
        )

        self.chunk_pl = self.chunk_pl.join(pvar_df, on="SNP", how="inner", suffix="_pvar")

        self.chunk_pl = self.chunk_pl.with_columns(
            pl.when(pl.col("REF") > pl.col("ALT"))
            .then(pl.col("BETA"))
            .otherwise(-pl.col("BETA")),
            pl.when(pl.col("REF") > pl.col("ALT"))
            .then(pl.col("AF"))
            .otherwise(1.0 - pl.col("AF"))
        )

    def ingest_data(self):
        """Append harmonized data to TileDB."""
        pl.Config.set_tbl_cols(-1) 
        print(self.chunk_pl)
        self.chunk_pl = self.chunk_pl.select(self.tiledb_types.keys())
        try:
            tiledb.from_pandas(
                uri=self.uri,
                dataframe=self.chunk_pl.to_pandas(),
                index_dims=self.dimension_tiledb,
                column_types=self.tiledb_types,
                mode="append",
            )
            logger.info(f"Successfully appended chunk to TileDB")
        except Exception as e:
            logger.error(f"Failed to append chunk to TileDB: {e}")
            raise

    def create_metadata(self, file_path: str, pvar_file: str = None):
        """Create and store metadata in TileDB."""
        metadata = {
            "file_path": file_path,
            "type_sumstat": self.type_sumstat,
        }

        if self.type_sumstat == "scqtl":
            celltype = self.chunk_pl["CELL"].unique().to_list()
            gene_list = self.chunk_pl["GENE"].unique().to_list()
            if not celltype:
                raise HarmonizationError("No cell types found in the data")
            metadata["cell_type"] = celltype
            metadata["genes"] = gene_list
        else:
            if "TRAIT" not in self.chunk_pl.columns:
                raise HarmonizationError("TRAIT column is missing in the data")
            if self.chunk_pl["TRAIT"].is_empty():
                raise HarmonizationError("TRAIT column is empty in the data")
            metadata["trait"] = self.chunk_pl["TRAIT"].unique().to_list()

        with tiledb.open(self.uri, mode='w') as array:
            array.meta["metadata"] = json.dumps(metadata)

        logger.info("Metadata created successfully")