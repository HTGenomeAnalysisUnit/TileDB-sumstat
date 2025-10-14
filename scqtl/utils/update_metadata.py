import tiledb
import pyarrow.compute as pc
import json
import pandas as pd
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

def update_metadata_single_cell_type(uri: str, cell: str) -> None:
    """
    Overwrite GENE_CELLTYPE metadata for a single cell type with unique genes from TileDB array.

    Args:
        uri: TileDB array URI.
        cell_type: Cell type to update (e.g., 'CD14__mono').
    """
    try:
        logger.info(f"Updating metadata for cell type {cell} in TileDB array at {uri}")

        # Open TileDB array
        with tiledb.open(uri, mode="r") as A:
            # Query GENE dimension for the specific cell type
            query = A.query(
                dims=["GENE"],
                attrs=[],
                cond=f"CELL == '{cell}'",
                return_arrow=True
            )
            # Get ChunkedArray for GENE
            gene_data = query.df[:, cell, :, :]["GENE"]

            # Compute unique genes
            unique_array = pc.unique(gene_data)
            unique_genes = unique_array.to_pandas().tolist()  # Convert to list
            logger.info(f"Found {len(unique_genes)} unique genes for {cell}")

            # Read existing metadata
            existing_gene_celltype = json.loads(A.meta.get("GENE_CELLTYPE", "{}"))

        # Update metadata for the cell type
        existing_gene_celltype[cell] = unique_genes
        logger.info(f"Updated metadata dictionary for {cell} with {len(unique_genes)} genes")

        # Write updated metadata back to TileDB
        with tiledb.open(uri, mode="w") as A:
            A.meta["GENE_CELLTYPE"] = json.dumps(existing_gene_celltype)
        logger.info(f"Successfully overwrote GENE_CELLTYPE metadata for {cell}")

    except Exception as e:
        logger.error(f"Failed to update metadata for {cell}: {e}")
        raise 