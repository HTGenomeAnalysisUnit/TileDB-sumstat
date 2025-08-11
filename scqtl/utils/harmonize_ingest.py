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
    def __init__(self, file_mapping: str, chunk_size: int, file_path: str, uri: str, pvar_file: str, type_sumstat: str):
        self.file_mapping = file_mapping
        self.chunk_size = chunk_size
        self.file_path = file_path
        self.uri = uri
        self.pvar_file = pvar_file
        self.type_sumstat = type_sumstat
        self.mapping_types = {}
        self.tiledb_types = {}
        self.dimension_tiledb = []
        self.chunk_pl = None  # will hold the harmonized dataframe

    def create_dtypes(self):
        """Create the mapping and dtype definitions."""
        df = pd.read_csv(self.file_mapping, sep=r"\s+", header=None, names=["key", "value"])
        self.mapping_types = dict(zip(df["key"], df["value"]))

        if self.type_sumstat == "gwas":
            self.tiledb_types = {
                "CHR": np.uint16,
                "TRAIT": str,
                "POS": np.uint32,
                "SNP": str,
                "DIST": np.int64,
                "AF": np.float32,
                "BETA": np.float32,
                "SE": np.float32,
                "P": np.float64,
                "N": np.int64,
            }
            self.dimension_tiledb = ["CHR", "TRAIT", "POS"]
        else:
            self.tiledb_types = {
                "CHR": np.uint16,
                "CELL": str,
                "GENE": str,
                "POS": np.uint32,
                "SNP": str,
                "DIST": np.int64,
                "AF": np.float32,
                "BETA": np.float32,
                "SE": np.float32,
                "P": np.float64,
                "N": np.int64,
            }
            self.dimension_tiledb = ["CHR", "CELL", "GENE", "POS"]

    def harmonize(self):
        """Load and rename columns, and ensure CHR/POS exist."""
        self.chunk_pl = pl.read_csv(
            self.file_path,
            separator="\t",
            columns=list(self.mapping_types.keys()),
            low_memory=True,
        )
        self.chunk_pl = self.chunk_pl.rename(self.mapping_types)
        self.chunk_pl = self.chunk_pl.with_columns(
            (pl.lit("chr") + pl.col("SNP"))
        )

        # If CHR/POS missing, extract from SNP ID
        if "CHR" not in self.chunk_pl.columns or "POS" not in self.chunk_pl.columns:
            self.chunk_pl = self.chunk_pl.with_columns(
                pl.col("SNP")
                .str.split_exact(":", 4)
                .struct.rename_fields(["CHR", "POS", "REF", "ALT"])
                .alias("fields")
            ).unnest("fields")

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

    def create_metadata(self):
        """Create and store metadata in TileDB."""
        metadata = {
            "file_path": self.file_path,
            "pvar_file": self.pvar_file,
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

    def run(self):
        """Run the full harmonization + ingestion pipeline."""
        self.create_dtypes()
        self.harmonize()
        self.align_alleles()
        self.ingest_data()
        self.create_metadata()