import tiledb
import logging
import json
import pandas as pd

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def _export_metadata_to_csv(
            uri: str = None,
            output_path: str = None):

    """Export metadata to CSV format."""
    if output_path is None:
        output_path = f"{uri}_metadata.csv"
    
        # Get the merged metadata from TileDB
        with tiledb.open(uri, "r") as array:
            merged_metadata_json = array.meta.get("merged_metadata", "{}")
    
        if not merged_metadata_json:
            logger.warning("No merged metadata found in TileDB array")
            return
    
        merged_metadata = json.loads(merged_metadata_json)
    
        # Create rows for CSV
        rows = []

        # Process QTL data (cell types)
        if "CELL" in merged_metadata and merged_metadata["CELL"]:
            for cell_type in merged_metadata["CELL"]:
                if cell_type in merged_metadata:
                    cell_data = merged_metadata[cell_type]
                    for chrom, genes in cell_data.items():
                        for gene_id, gene_metadata in genes.items():
                            row = {
                                "CHR": chrom,
                                "CELL": cell_type,
                                "GENE": gene_id,
                                "ACAT": gene_metadata.get("ACAT", ""),
                                "N": gene_metadata.get("N", ""),
                                "PHENO_VAR": gene_metadata.get("PHENO_VAR", ""),
                                "MIN_P": gene_metadata.get("MIN_P", "")
                                }
                            rows.append(row)
    
        # Process GWAS data (traits)
        if "traits" in merged_metadata and merged_metadata["traits"]:
            for trait in merged_metadata["traits"]:
                if trait in merged_metadata:
                    trait_data = merged_metadata[trait]
                    row = {
                    "TRAIT": trait,
                    "N": trait_data.get("N", ""),
                    "PHENO_VAR": trait_data.get("PHENO_VAR", ""),
                    "MIN_P": trait_data.get("MIN_P", ""),
                    "N_CASES": trait_data.get("N_CASES", ""),
                    "N_CONTROLS": trait_data.get("N_CONTROLS", "")
                    }
                    rows.append(row)
    
        # Create DataFrame and save to CSV
        if rows:
            df = pd.DataFrame(rows)
        
            # Reorder columns for better readability
            if "CELL" in df.columns:
                # QTL format
                column_order = ["CHR", "CELL", "GENE", "ACAT", "MIN_P", "N", "PHENO_VAR"]
                # Only include columns that exist in the DataFrame
                column_order = [col for col in column_order if col in df.columns]
                df = df[column_order]
            else:
                # GWAS format
                column_order = ["TRAIT", "N", "PHENO_VAR", "MIN_P", "N_CASES", "N_CONTROLS"]
                column_order = [col for col in column_order if col in df.columns]
                df = df[column_order]
        
            df.to_csv(output_path, index=False)
            logger.info(f"Metadata exported to {output_path}")
        
            # Print summary
            if "CELL" in df.columns:
                logger.info(f"Exported {len(df)} gene-cell type combinations")
                logger.info(f"Cell types: {df['CELL'].nunique()}")
                logger.info(f"Genes: {df['GENE'].nunique()}")
                logger.info(f"Chromosomes: {df['CHR'].nunique()}")
            else:
                logger.info(f"Exported {len(df)} traits")
        
            return df
        else:
            logger.warning("No metadata found to export")
            return pd.DataFrame()