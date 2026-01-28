import polars as pl
from pathlib import Path
import _align_alleles_with_pvar

def _align_alleles(
          chunk_pl:pl.DataFrame,
          pvar_file: Path
          ):
        #Remove the current SNPID
        if "SNPID" in chunk_pl.columns:
            chunk_pl = chunk_pl.drop("SNPID")
        swap = pl.col("A1") < pl.col("A2")
        #A2 is the effect allel while A is the non effect one        
        chunk_pl = chunk_pl.with_columns([
        pl.when(swap).then(pl.col("A1")).otherwise(pl.col("A2")).alias("EA"),
        # set NEA to the larger allele
        pl.when(swap).then(pl.col("A2")).otherwise(pl.col("A1")).alias("NEA")]
        )
        chunk_pl = chunk_pl.with_columns(
                pl.concat_str(
                    [
                    pl.col("CHR"),
                    pl.col("POS").cast(pl.Utf8),            # cast POS if numeric
                    pl.col("EA"),       # lexicographically smaller
                    pl.col("NEA")      # lexicographically larger
                    ],
                    separator=":"
                    ).alias("SNPID")
            )

        if pvar_file:
            _align_alleles_with_pvar(chunk_pl)
            chunk_pl.drop(["REF","ALT"])
        else:
            # flip the sign of BETA when swapping
            chunk_pl = chunk_pl.with_columns([
                pl.when(swap).then(-pl.col("BETA")).otherwise(pl.col("BETA")).alias("BETA"),
                # flip EAF to 1 - EAF when swapping
                pl.when(swap).then(1.0 - pl.col("EAF")).otherwise(pl.col("EAF")).alias("EAF"),
                # set EA to the smaller allele
            ])
        chunk_pl = chunk_pl.drop(["A1","A2"])
        chunk_pl = chunk_pl.with_columns(
                    pl.concat_str(
                        pl.lit("chr"),
                        pl.col("SNPID")).alias("SNPID"))
        return chunk_pl