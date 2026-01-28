import polars as pl
import tiledb
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def _ingest_data(
          tiledb_types: str = None,
          uri: str = None,
          dimension_tiledb: int = None,
          file_path: str = None,
          type_sumstat: str = None):
        """Append harmonized data to TileDB."""
        pl.Config.set_tbl_cols(-1)
        if type_sumstat == "gwas":
            dedup_keys = ["CHR", "POS", "TRAIT"]
        else:
            dedup_keys = ["CHR", "POS", "CELL", "GENE"]

        # Option A (recommended): window count -> keep only rows whose group count == 1
        chunk_pl = (
            chunk_pl
                .with_columns(pl.count().over(dedup_keys).alias("_grp_count"))
                .filter(pl.col("_grp_count") == 1)
                .drop("_grp_count")
                )
        chunk_pl = chunk_pl.with_columns([
            pl.col("CHR").cast(pl.UInt16),
            pl.col("POS").cast(pl.UInt32)
        ])
        chunk_pl_ingest = chunk_pl.select(tiledb_types.keys())
        chunk_pl_ingest = chunk_pl_ingest.drop_nulls()
        try:
            tiledb.from_pandas(
                uri=uri,
                dataframe=chunk_pl_ingest.to_pandas(),
                index_dims=dimension_tiledb,
                column_types=tiledb_types,
                allows_duplicates = False,
                mode="append"
            )
            logger.info(f"Successfully appended chunk to TileDB for file {file_path}")
        except Exception as e:
            logger.error(f"Failed to append chunk to TileDB for file {file_path}: {e}")
            raise