import os
from dask import delayed
import pyarrow as pa
import pyarrow.csv
import pyarrow.parquet as pq
import polars as pl
import numpy as np
import math
from scipy.stats import chi2, cauchy



@delayed
def batch_query_tiledb(array, queries, output_dir, batch_index):
    results = []
    for start, stop, cell_type in queries:
        result = array.query(return_arrow=True, dims=['gene', 'cell_type'], attrs=['p-value', 'beta', 'SNP']).df[start:stop, :, cell_type]
        results.append(result)
    
    # Concatenate all results into a single PyArrow Table
    combined_result = pa.concat_tables(results)
    
    # Check if the combined result has at least one row
    if combined_result.num_rows > 0:
        # Define the output file path
        output_path = os.path.join(output_dir, f"batch_{batch_index}.csv")
        
        # Write the combined result to a CSV file
        pyarrow.csv.write_csv(combined_result, output_path, write_options=pa.csv.WriteOptions(include_header=True))
        #gc.collect()
        return output_path
    else:
        # Return None if the combined result is empty
        #gc.collect()
        return None


def compute_pheno_variance(df):
    df_pl = pl.from_pandas(df)
    median_pl = df_pl.select(((pl.col("SE") ** 2) * pl.col("N") * 2 * pl.col("EAF") * (1 - pl.col("EAF"))).median())
    median_value = median_pl.item()
    median_str = str(median_value)
    return median_str

def acat_optimized(pvals_series: pl.Series, small: float = 1e-15) -> float:
    """
    Optimized ACAT implementation using vectorized operations with Polars.
    
    Parameters
    ----------
    pvals_series : pl.Series
        Series of p-values in (0, 1]. NaNs are ignored.
    small : float, optional
        Threshold below which we use the approximation.
        
    Returns
    -------
    float
        ACAT p-value in [0, 1].
    """
    # Convert to numpy array for vectorized operations
    p = pvals_series.to_numpy()
    
    # Remove NaNs
    valid_mask = ~np.isnan(p)
    if not np.any(valid_mask):
        return float("nan")
        
    p = p[valid_mask]
    
    # Validate p-values
    if np.any(p < 0) or np.any(p > 1):
        raise ValueError("All p-values must be in the interval [0, 1].")
    
    # Vectorized computation of transformed values
    t = np.empty_like(p)
    
    # Create masks
    small_mask = p < small
    regular_mask = ~small_mask
    
    # Apply transformations based on masks
    if np.any(regular_mask):
        t[regular_mask] = np.tan((0.5 - p[regular_mask]) * math.pi)
    
    if np.any(small_mask):
        t[small_mask] = 1.0 / (math.pi * p[small_mask])
    
    # Compute ACAT statistic
    cct_stat = np.mean(t)
    
    # Calculate ACAT p-value
    p_acat = cauchy.sf(cct_stat)
    
    # Clamp to [0,1] for numerical stability
    return max(0.0, min(p_acat, 1.0))


def z_to_p_via_chi2(z):
    """
    Convert z-score to p-value using chi-square distribution (df=1).

    Parameters:
        z (float): The z-score

    Returns:
        float: Two-tailed p-value
    """
    chi2_stat = z**2
    p_value = chi2.sf(chi2_stat, df=1)  # one-tailed
    return p_value
    