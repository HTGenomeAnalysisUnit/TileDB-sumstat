from tdbsumstat import HarmonizationError
import polars as pl
from tdbsumstat.utils import acat_optimized
from tdbsumstat.utils import compute_pheno_variance
import logging
import _write_individual_metadata

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def _create_metadata(
                    type_sumstat:str, 
                    type_trait:str,
                    file_path: str,
                    ):
        """Create and store metadata as individual JSON files."""
        metadata = {
            "traits": [],
            "CELL": []
        }
        chunk_pl = chunk_pl.drop_nulls()
        if type_sumstat == "qtl":
            # Get unique cell types
            celltypes = chunk_pl["CELL"].unique().to_list()
            if not celltypes:
                raise HarmonizationError("No cell types found in the data")

            # Update CELL list
            metadata["CELL"] = celltypes

            # Process each cell type
            for cell in celltypes:
                # Filter by this cell type and compute ACAT per gene
                df_cell = chunk_pl.filter(pl.col("CELL") == cell)
            
                # Group by CHR and GENE, compute ACAT for each group
                chr_gene_agg = df_cell.group_by(["CHR", "GENE"]).agg([
                    pl.col("P").map_batches(
                        lambda s: pl.Series([acat_optimized(s)]),
                        return_dtype=pl.Float64
                    ).alias("ACAT_LIST"),
                    pl.col("P").min().alias("min_P"),
                    pl.col("N").first().alias("N"),
                    pl.col("PHENO_VAR").first().alias("PHENO_VAR")
                ])
                chr_gene_agg = chr_gene_agg.with_columns(
                    pl.col("ACAT_LIST").list.first().alias("ACAT")
                )
            
                # Initialize cell structure if not exists
                if cell not in metadata:
                    metadata[cell] = {}
            
                # Populate metadata with chromosome -> gene structure
                for row in chr_gene_agg.iter_rows(named=True):
                    chrom = str(row["CHR"])  # Convert to string for consistency
                    gene = row["GENE"]
                    acat_val = row["ACAT"]
                    n_val = float(row["N"])
                    min_p = float(row["min_P"])
                    pheno_val = float(row["PHENO_VAR"])
                
                    if chrom not in metadata[cell]:
                        metadata[cell][chrom] = {}
                
                    gene_metadata = {
                        "ACAT": float(acat_val),
                        "N": n_val,
                        "PHENO_VAR": pheno_val,
                        "MIN_P":min_p
                    }
                
                    metadata[cell][chrom][gene] = gene_metadata
                    
        else:  # GWAS case
            if "TRAIT" not in chunk_pl.columns:
                raise HarmonizationError("TRAIT column is missing in the data")
        
            traits = chunk_pl["TRAIT"].unique().to_list()
            metadata["traits"] = traits
        
            for trait in traits:
                df_trait = chunk_pl.filter(pl.col("TRAIT") == trait)
                pheno_var = compute_pheno_variance(df_trait, type_trait)
                n_val = df_trait["N"].unique().to_list()[0]
                min_p = df_trait["P"].min().to_list()[0]
            
                trait_metadata = {
                    "N": float(n_val),
                    "PHENO_VAR": float(pheno_var),
                    "MIN_P": float(min_p)
                }
            
                if type_trait == "binary":
                    n_cases = df_trait["N_CASES"].unique().to_list()[0]
                    n_controls = df_trait["N_CONTROLS"].unique().to_list()[0]
                    trait_metadata.update({
                        "N_CASES": float(n_cases),
                        "N_CONTROLS": float(n_controls)
                    })
            
                metadata[trait] = trait_metadata

        # Write individual metadata file instead of updating TileDB
        _write_individual_metadata(metadata, file_path)
    
        logger.info(f"Individual metadata file created for {file_path}")