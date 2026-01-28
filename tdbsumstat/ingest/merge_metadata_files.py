import os
import logging
from pathlib import Path
import tiledb
import json
import _export_metadata_to_csv
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def _merge_metadata_files(
        uri: str = None
        ):
        """Merge all individual metadata files into final TileDB metadata."""
        metadata_dir = f"{uri}_metadata_parts"
    
        if not os.path.exists(metadata_dir):
            logger.warning(f"No metadata directory found at {metadata_dir}")
            return
    
        merged_metadata = {
        "traits": [],
        "CELL": []
        }
    
        # Process all individual metadata files
        metadata_files = list(Path(metadata_dir).glob("*.json"))
        logger.info(f"Found {len(metadata_files)} metadata files to merge")
    
        for metadata_file in metadata_files:
            try:
                with open(metadata_file, 'r') as f:
                    file_metadata = json.load(f)
            
                # Merge traits (for GWAS)
                if "traits" in file_metadata and file_metadata["traits"]:
                    current_traits = set(merged_metadata.get("traits", []))
                    new_traits = set(file_metadata["traits"])
                    merged_metadata["traits"] = list(current_traits.union(new_traits))
                
                    for trait in file_metadata["traits"]:
                        if trait in file_metadata:
                            if trait not in merged_metadata:
                                merged_metadata[trait] = file_metadata[trait]
                            else:
                                logger.warning(f"Trait {trait} already exists in metadata, overwriting")
                                merged_metadata[trait] = file_metadata[trait]
            
                # Merge CELL and cell metadata (for QTL)
                if "CELL" in file_metadata and file_metadata["CELL"]:
                    current_cells = set(merged_metadata.get("CELL", []))
                    new_cells = set(file_metadata["CELL"])
                    merged_metadata["CELL"] = list(current_cells.union(new_cells))
                
                    for cell in file_metadata["CELL"]:
                        if cell in file_metadata:
                            if cell not in merged_metadata:
                                merged_metadata[cell] = {}
                        
                            for chrom, genes in file_metadata[cell].items():
                                if chrom not in merged_metadata[cell]:
                                    merged_metadata[cell][chrom] = {}
                            
                                for gene, gene_data in genes.items():
                                    if gene in merged_metadata[cell][chrom]:
                                        logger.warning(f"Gene {gene} already exists in cell {cell} chromosome {chrom}, overwriting")
                                    merged_metadata[cell][chrom][gene] = gene_data
            
                logger.info(f"Processed {metadata_file.name}")
            
            except Exception as e:
                logger.error(f"Error processing metadata file {metadata_file}: {e}")
                continue
    
        # Store the final merged metadata in TileDB
        with tiledb.open(uri, "w") as array:
            array.meta["merged_metadata"] = json.dumps(merged_metadata)
    
        # Export to CSV
        _export_metadata_to_csv()
    
        # Log summary
        cell_count = len(merged_metadata.get("CELL", []))
        trait_count = len(merged_metadata.get("traits", []))
    
        total_genes = 0
        for cell_type in merged_metadata.get("CELL", []):
            if cell_type in merged_metadata:
                for chrom_data in merged_metadata[cell_type].values():
                    total_genes += len(chrom_data)
    
        logger.info(f"Final merged metadata: {cell_count} cell types, {trait_count} traits, {total_genes} total genes")
