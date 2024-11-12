import polars as pl
def process_write_chunk(chunk, SNP_list, file_stream):
    SNP_list_polars = pl.DataFrame(SNP_list)
    chunk_polars = pl.DataFrame(chunk)
    chunk_polars = chunk_polars.with_columns(
        pl.col("SNP").str.split_exact("_", 3)
        .struct.rename_fields(["pos",'A0','A1'])
        .alias("fields")
    ).unnest('fields')
    
    chunk_polars = chunk_polars.with_columns([
        pl.col("position").cast(pl.UInt32),
        pl.col("A0").cast(pl.String),
        pl.col("A1").cast(pl.String)
    ])

    # Perform the join operation with Polars
    subset_SNPs_merge = chunk_polars.join(
    SNP_list_polars,
    on=['position'],
    how="inner"
    )
    subset_SNPs_merge = subset_SNPs_merge.select('cell_type', 'gene', 'SNP', 'position', 'allele0', 'allele1', 'af', 'beta', 'se', 'p-value') 
    #Append the merged chunk to CSV
    subset_SNPs_merge.write_csv(file_stream)
