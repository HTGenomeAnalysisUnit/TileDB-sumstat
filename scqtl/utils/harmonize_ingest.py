import numpy as np
import polars as pl
import pandas as pd
import tiledb
import json
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class HarmonizationError(Exception):
    """Custom exception for harmonization errors."""
    pass

class Harmonize():
    def __init__(self, file_mapping: str, chunk_size: int, file_path: str, uri: str, pvar_file: str, type_sumstat: str):
        self.chunk_size = chunk_size
        self.file_path = file_path
        self.uri = uri
        self.pvar_file = pvar_file
        self.type_sumstat = type_sumstat
        self.file_mapping = file_mapping

    def create_dtypes(self):
        """Create the dtypes for the Polars DataFrame taking the fields from a file and converting them."""
        df = pd.read_csv("data.txt", sep="\s+", header=None, names=["key", "value"])
        # Convert to dictionary
        self.mapping_types = dict(zip(df["key"], df["value"]))  
        # Define the dtypes for TileDB
        # This is a dictionary mapping column names to their types
        # The types are defined based on the expected data types in the file
        # For example, "CHR" is an unsigned 16-bit integer, "POS" is an unsigned 32-bit integer, etc.
        # This mapping is used to ensure that the data is stored in the correct format in TileDB
        if self.type_sumstat =="gwas":
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
        chunk_pl = pl.read_csv(
            self.file_path,
            separator="\t",
            columns= self.mapping_types.keys(),
            low_memory=True,
        )
        # Change column name for each key value
        chunk_pl = chunk_pl.rename(self.mapping_types)
        chunk_pl = chunk_pl.with_columns(
                (pl.lit("chr") + pl.col("SNP"))               
            )
        # Create CHR POS REF ALT from variant ID if CHR POS are not present
        if "CHR" not in chunk_pl.columns or "POS" not in chunk_pl.columns:
            chunk_pl = chunk_pl.with_columns(
            pl.col("SNP")
            .str.split_exact(":", 4)
            .struct.rename_fields(["CHR", "POS","REF","ALT"])
            .alias("fields")
            ).unnest("fields")
        
        # Not sure about this
        #if type_sumstat == "scqtl":
        #    chunk_pl = chunk_pl.with_columns(
        #        pl.when(pl.col("start_distance").str.contains("vs"))
        #        .then(None)
        #        .otherwise(pl.col("start_distance"))
        #        .cast(pl.Int64, strict=False)
        #        .alias("DIST"),
        #    )

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
                new_columns=["CHROM", "POS", "SNP", "REF", "ALT", "INFO"],  # Rename columns for clarity
                dtypes={"CHROM": pl.Utf8, "POS": pl.Utf8, "SNP": pl.Utf8, "REF": pl.Utf8, "ALT": pl.Utf8, "INFO":pl.Utf8},
            )
            chunk_pl = chunk_pl.join(pvar_df, on="SNP", how="inner", suffix="_pvar")        
            chunk_pl = chunk_pl.with_columns(
                pl.when(pl.col("REF") > pl.col("ALT"))
                .then(pl.col("BETA"))  # Keep BETA as is
                .otherwise(-pl.col("BETA")),  # Invert BETA
                pl.when(pl.col("REF") > pl.col("ALT"))
                .then(pl.col("AF"))  # Keep AF as is
                .otherwise(1.0 - pl.col("AF"))  # Compute 1 - AF
            )
    
        def ingest_data(self):
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
            """Create metadata for the TileDB array."""
            metadata = {
                "file_path": self.file_path,
                "pvar_file": self.pvar_file,
                "type_sumstat": self.type_sumstat,
            }
            if self.type_sumstat == "scqtl":
                celltype = chunk_pl["CELL"].unique().to_list()
                gene_list = chunk_pl["GENE"].unique().to_list()
                if not celltype:
                    raise HarmonizationError("No cell types found in the data")
                if metadata.get("cell_type") is None:
                    # Initialize cell_type if it doesn't exist
                    metadata["cell_type"] = celltype
                else:
                    # Append cell types to existing list
                    metadata["cell_type"] += celltype
                if not gene_list:
                    raise HarmonizationError("No genes found in the data")
                if metadata.get(celltype) is None:
                    # Initialize gene if it doesn't exist
                    metadata[celltype] = gene_list
                else:
                    # Append genes to existing list
                    metadata[celltype] += gene_list
            else:
                # For GWAS, we just need the TRAIT
                if "TRAIT" not in chunk_pl.columns:
                    raise HarmonizationError("TRAIT column is missing in the data")
                if chunk_pl["TRAIT"].is_empty():
                    raise HarmonizationError("TRAIT column is empty in the data")
                if metadata.get("trait") is None:
                    # Initialize trait if it doesn't exist
                    metadata["trait"] = chunk_pl["TRAIT"].unique().to_list()
                else:
                    # Append traits to existing list
                    metadata["trait"] += chunk_pl["TRAIT"].unique().to_list()
            
            # Write metadata to TileDB array
            with tiledb.open(self.uri, mode='w') as array:
                array.meta["metadata"] = json.dumps(metadata)
            logger.info("Metadata created successfully")