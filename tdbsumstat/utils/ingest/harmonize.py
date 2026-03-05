"""Data harmonization – column renaming, allele standardisation, p-value computation."""
import logging
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats

logger = logging.getLogger(__name__)


class HarmonizeMixin:
    """Mixin that normalises a summary statistics Polars DataFrame.

    After calling :meth:`harmonize`, the processed data is stored in
    ``self.chunk_pl`` and column-type metadata is stored in
    ``self.tiledb_types``.
    """

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def harmonize(
        self,
        sumstat: pl.DataFrame,
        trait: str = None,
        cell: str = None,
        gene: str = None,
        pheno_var: float = None,
        n: int = None,
        n_controls: int = None,
        n_cases: int = None,
    ) -> None:
        """Rename columns, compute missing fields, and set TileDB types.

        Parameters
        ----------
        sumstat:
            Raw Polars DataFrame to harmonise.
        trait:
            Trait identifier (GWAS only).
        cell:
            Cell-type identifier (QTL only).
        gene:
            Gene identifier (QTL only).
        pheno_var:
            Phenotypic variance (quant trait).
        n:
            Sample size.
        n_controls:
            Number of controls (binary GWAS).
        n_cases:
            Number of cases (binary GWAS).

        Raises
        ------
        HarmonizationError
            When required columns are missing and cannot be inferred.
        """
        from tdbsumstat.utils.ingest.errors import HarmonizationError

        self.chunk_pl = sumstat.rename(self.mapping_types)

        # Extract CHR/POS from SNPID if not present
        if "CHR" not in self.chunk_pl.columns or "POS" not in self.chunk_pl.columns:
            self.chunk_pl = self.chunk_pl.with_columns(
                pl.col("SNPID")
                .str.split_exact(":", 4)
                .struct.rename_fields(["CHR", "POS", "A1", "A2"])
                .alias("fields")
            ).unnest("fields")

        if "SNPID" in self.chunk_pl.columns:
            self.chunk_pl = self.chunk_pl.drop("SNPID")
        if "EAF" not in self.chunk_pl.columns:
            self.chunk_pl = self.chunk_pl.with_columns(pl.lit(0).alias("EAF"))
        if "DIST" not in self.chunk_pl.columns:
            self.chunk_pl = self.chunk_pl.with_columns(pl.lit(1).alias("DIST"))

        # Optional MAF filter
        if self.maf is not None:
            self.chunk_pl = self.chunk_pl.with_columns(
                pl.min_horizontal(pl.col("EAF"), 1 - pl.col("EAF")).alias("MAF")
            ).filter(pl.col("MAF") >= self.maf)

        # Sample-size and phenotypic-variance handling
        self._apply_sample_size(n, n_cases, n_controls, pheno_var)

        # Second SNPID-drop guard (in case rename created one)
        if "SNPID" in self.chunk_pl.columns:
            self.chunk_pl = self.chunk_pl.drop("SNPID")

        # Allele standardisation: A1 < A2 order
        self._standardise_alleles()

        # Assign RSID placeholder if absent
        if "RSID" not in self.chunk_pl.columns:
            self.chunk_pl = self.chunk_pl.with_columns(pl.lit("None").alias("RSID"))

        # Convert LOG10P → P
        if "LOG10P" in self.chunk_pl.columns:
            self.chunk_pl = self.chunk_pl.with_columns(
                (10 ** (-pl.col("LOG10P"))).alias("P")
            )

        # Compute / recalculate p-value
        if self.permuted:
            self.chunk_pl = self.chunk_pl.with_columns(
                (
                    pl.col("BETA").pow(2)
                    / pl.col("P").map_batches(
                        lambda x: pl.Series(stats.chi2.isf(x.to_numpy(), df=1))
                    )
                ).sqrt().alias("SE")
            )
        else:
            self.chunk_pl = self.chunk_pl.drop("P")
            self.chunk_pl = self.chunk_pl.with_columns(
                (pl.col("BETA") / pl.col("SE"))
                .pow(2)
                .map_batches(
                    lambda x: pl.Series(stats.chi2.sf(x.to_numpy(), df=1)),
                    return_dtype=pl.Float64,
                )
                .alias("P")
            )

        # Set dimension labels and TileDB type maps
        self._set_tiledb_types(trait=trait, cell=cell, gene=gene)

    def align_alleles(self) -> None:
        """Align effect/other alleles using an external pvar reference file.

        Requires ``self.pvar_file`` to be set and point to a valid file.

        Raises
        ------
        HarmonizationError
            If ``pvar_file`` is not provided.
        FileNotFoundError
            If the pvar file path does not exist.
        """
        from tdbsumstat.utils.ingest.errors import HarmonizationError

        if not self.pvar_file:
            raise HarmonizationError("pvar_file must be provided to verify alleles order")
        if not Path(self.pvar_file).is_file():
            raise FileNotFoundError(f"pvar_file {self.pvar_file} does not exist")

        pvar_df = pl.read_csv(
            self.pvar_file,
            separator="\t",
            has_header=True,
            schema_overrides={"CHROM": pl.Utf8, "POS": pl.Utf8, "SNPID": pl.Utf8,
                              "REF": pl.Utf8, "ALT": pl.Utf8},
        )
        self.chunk_pl = self.chunk_pl.join(pvar_df, on="SNPID", how="inner", suffix="_pvar")

        swap = pl.col("ALT") < pl.col("REF")
        self.chunk_pl = self.chunk_pl.with_columns([
            pl.when(swap).then(pl.col("BETA")).otherwise(-pl.col("BETA")),
            pl.when(swap).then(pl.col("EAF")).otherwise(1.0 - pl.col("EAF")),
        ])

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _apply_sample_size(self, n, n_cases, n_controls, pheno_var):
        """Add N, N_CASES, N_CONTROLS columns and apply MAC filter."""
        from tdbsumstat.utils.ingest.errors import HarmonizationError

        if self.type_trait == "quant":
            if "N" not in self.chunk_pl.columns:
                if n is not None:
                    self.chunk_pl = self.chunk_pl.with_columns(pl.lit(int(n)).alias("N"))
                else:
                    raise HarmonizationError("N column is missing and N parameter is not provided")
            if self.mac is not None:
                self.chunk_pl = self.chunk_pl.with_columns(
                    (2 * pl.col("N") * pl.min_horizontal(pl.col("EAF"), 1 - pl.col("EAF"))).alias("MAC")
                ).filter(pl.col("MAC") >= self.mac)
            self.chunk_pl = self.chunk_pl.with_columns(pl.lit(pheno_var).alias("PHENO_VAR"))

        elif self.type_trait == "binary":
            if not all(col in self.chunk_pl.columns for col in ["N_CASES", "N_CONTROLS"]):
                if n_cases is not None and n_controls is not None:
                    self.chunk_pl = self.chunk_pl.with_columns(
                        pl.lit(float(n_cases)).alias("N_CASES"),
                        pl.lit(float(n_controls)).alias("N_CONTROLS"),
                        pl.lit(float(n_cases) + float(n_controls)).alias("N"),
                    )
                else:
                    raise HarmonizationError("n_cases and n_controls columns are missing and were not provided")
                if self.mac is not None:
                    self.chunk_pl = self.chunk_pl.with_columns(
                        (2 * pl.col("N") * pl.min_horizontal(pl.col("EAF"), 1 - pl.col("EAF"))).alias("MAC")
                    ).filter(pl.col("MAC") >= self.mac)
        else:
            from tdbsumstat.utils.ingest.errors import HarmonizationError
            raise HarmonizationError("Type of trait must be either binary or quant")

    def _standardise_alleles(self):
        """Sort A1/A2 lexicographically; flip BETA and EAF when swapped."""
        swap = pl.col("A1") < pl.col("A2")

        self.chunk_pl = self.chunk_pl.with_columns([
            pl.when(swap).then(pl.col("A1")).otherwise(pl.col("A2")).alias("EA"),
            pl.when(swap).then(pl.col("A2")).otherwise(pl.col("A1")).alias("NEA"),
        ])

        self.chunk_pl = self.chunk_pl.with_columns(
            pl.concat_str(
                [
                    pl.col("CHR"),
                    pl.col("POS").cast(pl.Utf8),
                    pl.col("EA"),
                    pl.col("NEA"),
                ],
                separator=":",
            ).alias("SNPID")
        )

        if self.pvar_file:
            self.align_alleles()
            self.chunk_pl = self.chunk_pl.drop(["REF", "ALT"])
        else:
            self.chunk_pl = self.chunk_pl.with_columns([
                pl.when(swap).then(-pl.col("BETA")).otherwise(pl.col("BETA")).alias("BETA"),
                pl.when(swap).then(1.0 - pl.col("EAF")).otherwise(pl.col("EAF")).alias("EAF"),
            ])

        self.chunk_pl = self.chunk_pl.drop(["A1", "A2"])
        self.chunk_pl = self.chunk_pl.with_columns(
            pl.concat_str(pl.lit("chr"), pl.col("SNPID")).alias("SNPID")
        )

    def _set_tiledb_types(self, trait, cell, gene):
        """Populate ``self.tiledb_types`` and add dimension columns."""
        if self.type_sumstat == "gwas":
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
                self.chunk_pl = self.chunk_pl.with_columns(pl.lit(trait).alias("TRAIT"))
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
                self.chunk_pl = self.chunk_pl.with_columns(pl.lit(cell).alias("CELL"))
            if "GENE" not in self.chunk_pl.columns:
                self.chunk_pl = self.chunk_pl.with_columns(pl.lit(gene).alias("GENE"))
