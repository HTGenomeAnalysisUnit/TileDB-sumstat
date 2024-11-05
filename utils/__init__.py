import os
from dask import delayed, compute
import pyarrow as pa
import pyarrow.csv
import pyarrow.parquet as pq
import gc
import scipy.stats as stats

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