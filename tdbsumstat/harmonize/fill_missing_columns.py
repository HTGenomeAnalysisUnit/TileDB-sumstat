import polars as pl
from tdbsumstat import HarmonizationError
import numpy as np
import _create_mapping
from pathlib import Path


def _fill_missing_columns(
        chunk_pl: pl.DataFrame,
        mapping_file:Path,
        type_trait: str,
        type_sumstat: str,
        mac: int | None = 10, 
        trait: str | None = None, 
        cell: str | None = None, 
        gene: str | None = None, 
        pheno_var:int | None = None, 
        n: int | None = None, 
        n_controls: int | None = None,
        n_cases: int | None = None):

        """Load and rename columns, and ensure CHR/POS exist."""
        chunk_pl = _create_mapping(sumstat = chunk_pl, mapping_file = mapping_file)
        if "CHR" not in chunk_pl.columns or "POS" not in chunk_pl.columns:
            chunk_pl = chunk_pl.with_columns(
                pl.col("SNPID")
                .str.split_exact(":", 4)
                .struct.rename_fields(["CHR", "POS", "A1", "A2"])
                .alias("fields")
            ).unnest("fields")
    
        if "EAF" not in chunk_pl.columns:
            chunk_pl = chunk_pl.with_columns(pl.lit(0).alias("EAF"))
        if "DIST" not in chunk_pl.columns:
            chunk_pl = chunk_pl.with_columns(pl.lit(1).alias("DIST"))
        #Check for removing double headers

        if type_trait == "quant": 
            if not "N" in chunk_pl.columns:
                if n is not None:
                    chunk_pl = chunk_pl.with_columns(
                    pl.lit(int(n)).alias("N")
                    )
                else:
                    raise HarmonizationError("N column is missing and N parameter is not provided")
            if mac is not None:
                chunk_pl = chunk_pl.with_columns(
                (2 * pl.col("N") * pl.min_horizontal(pl.col("EAF"), 1 - pl.col("EAF")))
                .alias("MAC")
                 ).filter(pl.col("MAC") >= mac)
            
            chunk_pl = chunk_pl.with_columns(
                    pl.lit(pheno_var).alias("PHENO_VAR")
                )
        elif type_trait == "binary": 
            if not all(sample_size in chunk_pl.columns for sample_size in ["N_CASES", "N_CONTROLS"]):
                if not None in [n_cases, n_controls]:
                    chunk_pl = chunk_pl.with_columns(
                    pl.lit(float(n_cases)).alias("N_CASES"),
                    pl.lit(float(n_controls)).alias("N_CONTROLS"),
                    pl.lit(float(n_cases) + float(n_controls)).alias("N"),
                    )
                else:
                    raise HarmonizationError("n_cases and n_controls columns are missing and were not provided")
                if mac is not None:
                    chunk_pl = chunk_pl.with_columns(
                    (2 * pl.col("N") * pl.min_horizontal(pl.col("EAF"), 1 - pl.col("EAF")))
                    .alias("MAC")
                    ).filter(pl.col("MAC") >= mac)
        else:
            raise HarmonizationError("Type of trait must be either binary or quant")
                
        if type_sumstat=="gwas":
            tiledb_types = {
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
            if "TRAIT" not in chunk_pl.columns:
                chunk_pl = chunk_pl.with_columns(
                    pl.lit(trait).alias("TRAIT")
                )
        else:
            tiledb_types = {
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
            if "CELL" not in chunk_pl.columns:
                chunk_pl = chunk_pl.with_columns(
                    pl.lit(cell).alias("CELL")
                )
            if "GENE" not in chunk_pl.columns:
                chunk_pl = chunk_pl.with_columns(
                    pl.lit(gene).alias("GENE")
                )
        if "RSID" not in chunk_pl.columns:
                chunk_pl = chunk_pl.with_columns(
                    pl.lit("None").alias("RSID")
                )
        if "LOG10P" in chunk_pl.columns:
            chunk_pl = chunk_pl.with_columns(
                    (10 ** (-pl.col("LOG10P"))).alias("P")
                )        