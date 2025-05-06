import numpy as np
import polars as pl
import tiledb
import json
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

def harmonize_ingest_data(chunk_size: int, file: str, uri: str, pvar_file: str) -> None:
    """
    Process genetic data in chunks, store in TileDB, and update metadata once.

    Args:
        chunk_size: Number of rows per chunk.
        file: Comma-separated string with file path and cell type (e.g., "path.tsv,cell_type").
        uri: TileDB array URI.
    """
    try:
        # Validate input
        file_parts = file.split(",")
        if len(file_parts) != 2:
            raise ValueError("Expected file input as 'filepath,cell_type'")
        file_path, cell_type = file_parts
        if not Path(file_path).is_file():
            raise FileNotFoundError(f"File {file_path} does not exist")

        logger.info(f"Processing file {file_path} with cell type {cell_type}")

        # Define schema
        dtypes = {
            "variant_id": pl.Utf8,
            "start_distance": pl.Utf8,
            "phenotype_id": pl.Utf8,
            "slope": pl.Float32,
            "slope_se": pl.Float32,
            "af": pl.Float32,
            "pval_nominal": pl.Float64,
            "ma_count": pl.Int64,
        }
        columns = list(dtypes.keys())

        # TileDB column types
        tiledb_types = {
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
        
        chunk_pl = pl.read_csv(
            file_path,
            separator="\t",
            columns=columns,
            dtypes=dtypes,
            low_memory=True,
        )
            # Clean start_distance
        chunk_pl = chunk_pl.with_columns(
            pl.col("variant_id")
            .str.split_exact(":", 4)
            .struct.rename_fields(["CHR", "POS","REF","ALT"])
            .alias("fields")
        ).unnest("fields").drop(["REF", "ALT"])

        chunk_pl = chunk_pl.with_columns(
                pl.when(pl.col("start_distance").str.contains("vs"))
                .then(None)
                .otherwise(pl.col("start_distance"))
                .cast(pl.Int64, strict=False)
                .alias("DIST"),
                pl.lit(None, pl.Utf8).alias("RSID"),
            )        
        # Create SNP identifier
        chunk_pl = chunk_pl.with_columns(
                pl.col("CHR").str.replace(r"chr", "").cast(pl.Utf8),
                pl.col("POS").cast(pl.Utf8),
                pl.when(pl.col("af") > 0)
                .then(pl.col("ma_count") / (pl.col("af") * 2))
                .otherwise(None)
                .cast(pl.Int64)
                .alias("N"),
            )
        
        pvar_df = pl.read_csv(
            pvar_file,
            separator="\t",
            has_header=True,
            new_columns=["CHROM", "POS", "SNP", "REF", "ALT", "INFO"],  # Rename columns for clarity
            dtypes={"CHROM": pl.Utf8, "POS": pl.Utf8, "SNP": pl.Utf8, "REF": pl.Utf8, "ALT": pl.Utf8, "INFO":pl.Utf8},
        )
        
        chunk_pl = chunk_pl.join(pvar_df, right_on="SNP", left_on="variant_id", how="inner", suffix="_pvar")        
        chunk_pl = chunk_pl.with_columns(
            pl.when(pl.col("REF") > pl.col("ALT"))
            .then(pl.col("slope"))  # Keep BETA as is
            .otherwise(-pl.col("slope"))  # Invert BETA
            .alias("BETA"),
            pl.when(pl.col("REF") > pl.col("ALT"))
            .then(pl.col("af"))  # Keep AF as is
            .otherwise(1.0 - pl.col("af"))  # Compute 1-AF
            .alias("AF"),
        )
        
        chunk_pl = chunk_pl.with_columns(
                (pl.lit("chr") + pl.col("variant_id")).alias("SNP"),
                pl.col("phenotype_id").alias("GENE"),
                pl.col("slope_se").alias("SE"),
                pl.col("pval_nominal").alias("P"),
                pl.lit(cell_type).alias("CELL"),
            )

            # Select final columns
        chunk_processed = chunk_pl.select([
                "CHR", "CELL", "GENE", "POS", "SNP", "RSID", "DIST", "AF", "BETA", "SE", "P", "N"
            ])

            # Convert to pandas for TileDB
        chunk_df = chunk_processed.to_pandas()

            # Validate DataFrame
        if not all(col in chunk_df.columns for col in tiledb_types):
                raise ValueError(f"Chunk missing required columns")

            # Write to TileDB
        try:
                tiledb.from_pandas(
                    uri=uri,
                    dataframe=chunk_df,
                    index_dims=["CHR", "CELL", "GENE", "POS"],
                    column_types=tiledb_types,
                    mode="append",
                )
                logger.info(f"Successfully appended chunk to TileDB")
        except Exception as e:
                logger.error(f"Failed to append chunk to TileDB: {e}")
                raise

    except Exception as e:
        logger.error(f"Error processing data: {e}")
        raise