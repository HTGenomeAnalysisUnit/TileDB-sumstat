from .exceptions import HarmonizationError
from pathlib import Path
import polars as pl

def _align_alleles_with_pvar(
          pvar_file: Path | None = None):
        """Align alleles based on pvar file."""
        if pvar_file is None:
            raise HarmonizationError("pvar_file must be provided to verify alleles order")
        if not pvar_file.is_file():
            raise FileNotFoundError(f"pvar_file {pvar_file} does not exist")

        #The ALT must correspond to the alternative allele
        pvar_df = pl.read_csv(
            pvar_file,
            separator="\t",
            has_header=True,
            dtypes={"CHROM": pl.Utf8, "POS": pl.Utf8, "SNPID": pl.Utf8,
                    "REF": pl.Utf8, "ALT": pl.Utf8},
        )
        #Here we assume the SNPID is alphabetically sortedin both pvar and summary statistics
        chunk_pl = chunk_pl.join(pvar_df, on="SNPID", how="inner", suffix="_pvar")

        swap = pl.col("ALT")< pl.col("REF")

        chunk_pl = chunk_pl.with_columns([
            pl.when(swap).then(pl.col("BETA")).otherwise(-pl.col("BETA")),
            pl.when(swap).then(pl.col("EAF")).otherwise(1.0-pl.col("EAF"))
            ])
