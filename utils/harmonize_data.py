import numpy as np
import pandas as pd
import scipy.stats as stats

def calculate_p_values(beta, se):
    """Calculate two-tailed p-values from t-statistics."""
    return stats.t.sf(np.abs(beta / se), df=se.size - 1) * 2

def harmonize_data(chunk, cell_type):
    # Split the SNP into chrom, pos, ref, alt directly without row-wise apply
    chrompos_split = chunk['variant_id'].str.split("_", expand=True)
    chrom = chrompos_split[0]
    pos = chrompos_split[1]
    ref = chrompos_split[2]
    alt = chrompos_split[3]

    # Use vectorized numpy operations for fast processing
    allele0 = np.where(ref < alt, ref, alt)
    allele1 = np.where(ref >= alt, ref, alt)

    # Create sorted SNP and new beta values
    sorted_snp = chrom + "_" + pos + "_" + allele0 + "_" + allele1
    new_beta = np.where(allele0 == alt, -chunk['slope'], chunk['slope'])
    new_af = np.where(allele0 == alt, 1-chunk['af'], chunk['af'])

    # Assign results back to the chunk
    chunk['SNP'] = sorted_snp
    chunk['gene'] = chunk["phenotype_id"]
    chunk['position'] = pos
    chunk['allele1'] = allele0
    chunk['allele0'] = allele1
    chunk['beta'] = new_beta
    chunk['se'] = chunk['slope_se']
    chunk['cell_type'] = cell_type
    chunk['af'] = new_af
    chunk['p-value'] = calculate_p_values(chunk['beta'],chunk['beta'])
    # Select necessary columns for TileDB
    chunk_processed = chunk[['SNP', 'gene', 'position', 'allele0', 'allele1', 'af', 'beta', 'se', 'p-value', 'cell_type']]
    return chunk_processed

