"""TileDB data writer – appends harmonised data to an existing TileDB array."""
import logging

import polars as pl
import tiledb

logger = logging.getLogger(__name__)


class WriterMixin:
    """Mixin that deduplicates and writes ``self.chunk_pl`` into TileDB."""

    def ingest_data(self, file_path: str) -> None:
        """Deduplicate and append ``self.chunk_pl`` to the TileDB array.

        Rows with duplicate index-dimension keys are removed before writing.
        Only the columns declared in ``self.tiledb_types`` are written.

        Parameters
        ----------
        file_path:
            Original file path; used only for logging messages.

        Raises
        ------
        Exception
            Re-raises any exception from ``tiledb.from_pandas``.
        """
        pl.Config.set_tbl_cols(-1)

        dedup_keys = ["CHR", "POS", "TRAIT"] if self.type_sumstat == "gwas" else ["CHR", "POS", "CELL", "GENE"]

        # Drop ALL rows that share an index key with at least one other row.
        # Rows with ambiguous (duplicated) positions are treated as unreliable
        # and excluded from the TileDB write entirely.
        self.chunk_pl = (
            self.chunk_pl
            .with_columns(pl.len().over(dedup_keys).alias("_grp_count"))
            .filter(pl.col("_grp_count") == 1)
            .drop("_grp_count")
        )
        self.chunk_pl = self.chunk_pl.with_columns([
            pl.col("CHR").cast(pl.UInt16),
            pl.col("POS").cast(pl.UInt32),
        ])

        chunk_pl_ingest = self.chunk_pl.select(self.tiledb_types.keys()).drop_nulls()

        try:
            tiledb.from_pandas(
                uri=self.uri,
                dataframe=chunk_pl_ingest.to_pandas(),
                index_dims=self.dimension_tiledb,
                column_types=self.tiledb_types,
                allows_duplicates=False,
                mode="append",
            )
            logger.info("Successfully appended chunk to TileDB for file %s", file_path)
        except Exception as exc:
            logger.error("Failed to append chunk to TileDB for file %s: %s", file_path, exc)
            raise
