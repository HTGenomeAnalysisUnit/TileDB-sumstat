from scipy import stats
import polars as pl
def _compute_pval():
    chunk_pl = chunk_pl.with_columns(
    (pl.col("BETA") / pl.col("SE")).pow(2).map_batches(
    lambda x: pl.Series(stats.chi2.sf(x.to_numpy(), df=1)),
    return_dtype=pl.Float64
    ).alias('P')
    )