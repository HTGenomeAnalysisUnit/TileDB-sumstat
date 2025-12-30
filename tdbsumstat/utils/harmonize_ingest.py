import logging
import json
from pathlib import Path
import os
import pandas as pd
import polars as pl
import numpy as np
from scipy import stats
import tiledb
import gwaslab as gl
from collections import defaultdict
from tdbsumstat.utils import acat_optimized, compute_pheno_variance


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


class HarmonizationError(Exception):
    """Custom exception for harmonization errors."""
    pass


class Harmonize:
    def __init__(self, mapping_file: str, uri: str, type_sumstat: str, pvar_file: str, type_trait: str, mac: int):
        self.mapping_file = mapping_file
        self.uri = uri
        self.pvar_file = pvar_file
        self.type_sumstat = type_sumstat
        self.type_trait = type_trait
        self.mapping_types = {}
        self.tiledb_types = {}
        self.dimension_tiledb = []
        self.mac = mac
    
    def create_mapping(self):
        df = pd.read_csv(self.mapping_file, header=None, names=["key", "value"])
        if df.empty:
            raise HarmonizationError("Mapping file is empty or not formatted correctly.")
        self.mapping_types = dict(zip(df["key"], df["value"]))
        #check that "BETA", "SE" are in the vlaues of the mapping_types
        if not all(col in self.mapping_types.values() for col in ["BETA", "SE"]):
            raise HarmonizationError("Mapping file must contain BETA and SE columns.")
        # Check if CHR and POS or SNP are present
        if not all(col for col in ["CHR", "POS"] if col in self.mapping_types.values()):
            if "SNPID" not in self.mapping_types.values():
                raise HarmonizationError("Mapping file must contain either CHR, and POS or SNPID columns.")
        
    def create_tiledb(self):
        """Create the mapping and dtype definitions."""
        pos_domain = (1, 300000000)  # Example range for genomic positions
        chr_domain = (1, 24)  # Example range for genomic positions
        attrs=[
                tiledb.Attr(name="SNPID", dtype="ascii", filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="RSID", dtype="ascii", filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="EAF", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="BETA", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="SE", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="P", dtype=np.float64, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
            ]
    
        if self.type_sumstat == "gwas":
            
            self.dimension_tiledb = ["CHR", "TRAIT", "POS"]
            dom = tiledb.Domain(
            tiledb.Dim(name="CHR", domain = chr_domain, dtype=np.uint16,  filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            tiledb.Dim(name="TRAIT", dtype="ascii",  var=True, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]) ),
            tiledb.Dim(name="POS", domain = pos_domain, dtype=np.uint32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
            )
        else:
            
            self.dimension_tiledb = ["CHR", "CELL", "GENE", "POS"]
            dom = tiledb.Domain(
                tiledb.Dim(name="CHR", domain = chr_domain, dtype=np.uint16,  filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Dim(name="CELL", dtype="ascii",  var=True, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]) ),
                tiledb.Dim(name="GENE",dtype="ascii", var=True, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]) ),
                tiledb.Dim(name="POS", domain = pos_domain, dtype=np.uint32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
                )
            
            attrs = attrs + [
            tiledb.Attr(name="DIST", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
            ]
        schema = tiledb.ArraySchema(
                domain=dom,
                attrs=attrs,
                sparse=True,
                allows_duplicates=False
        )
        tiledb.Array.create(self.uri, schema)
    
    def align_alleles(self):
        """Align alleles based on pvar file."""
        if not self.pvar_file:
            raise HarmonizationError("pvar_file must be provided to verify alleles order")
        if not Path(self.pvar_file).is_file():
            raise FileNotFoundError(f"pvar_file {self.pvar_file} does not exist")

        #The ALT must correspond to the alternative allele
        pvar_df = pl.read_csv(
            self.pvar_file,
            separator="\t",
            has_header=True,
            dtypes={"CHROM": pl.Utf8, "POS": pl.Utf8, "SNPID": pl.Utf8,
                    "REF": pl.Utf8, "ALT": pl.Utf8},
        )
        #Here we assume the SNPID is alphabetically sortedin both pvar and summary statistics
        self.chunk_pl = self.chunk_pl.join(pvar_df, on="SNPID", how="inner", suffix="_pvar")

        swap = pl.col("ALT")< pl.col("REF")

        self.chunk_pl = self.chunk_pl.with_columns([
            pl.when(swap).then(pl.col("BETA")).otherwise(-pl.col("BETA")),
            pl.when(swap).then(pl.col("EAF")).otherwise(1.0-pl.col("EAF"))
            ])

    def harmonize(self, 
                  sumstat, 
                  mac: int, 
                  trait: str = None, 
                  cell: str = None, 
                  gene: str = None, 
                  pheno_var:int = None, 
                  n: int = None, 
                  n_controls: int = None,
                  n_cases: int = None):
        """Load and rename columns, and ensure CHR/POS exist."""
        self.chunk_pl = sumstat.rename(self.mapping_types)
        # If CHR/POS missing, extract from SNP ID
        if "CHR" not in self.chunk_pl.columns or "POS" not in self.chunk_pl.columns:
            self.chunk_pl = self.chunk_pl.with_columns(
                pl.col("SNPID")
                .str.split_exact(":", 4)
                .struct.rename_fields(["CHR", "POS", "A1", "A2"])
                .alias("fields")
            ).unnest("fields")
    
        if "SNPID" in self.chunk_pl.columns:
            self.chunk_pl = self.chunk_pl.drop("SNPID")
        if "EAF" not in self.chunk_pl.columns:
            self.chunk_pl = self.chunk_pl.with_columns(pl.lit(0).alias("EAF"))
        if "DIST" not in self.chunk_pl.columns:
            self.chunk_pl = self.chunk_pl.with_columns(pl.lit(1).alias("DIST"))
        #Check for removing double headers

        if self.type_trait == "quant": 
            if not "N" in self.chunk_pl.columns:
                if n is not None:
                    self.chunk_pl = self.chunk_pl.with_columns(
                    pl.lit(int(n)).alias("N")
                    )
                else:
                    raise HarmonizationError("N column is missing and N parameter is not provided")
            if self.mac is not None:
                self.chunk_pl = self.chunk_pl.with_columns(
                (2 * pl.col("N") * pl.min_horizontal(pl.col("EAF"), 1 - pl.col("EAF")))
                .alias("MAC")
                 ).filter(pl.col("MAC") >= self.mac)
            
            self.chunk_pl = self.chunk_pl.with_columns(
                    pl.lit(pheno_var).alias("PHENO_VAR")
                )
        elif self.type_trait == "binary": 
            if not all(sample_size in self.chunk_pl.columns for sample_size in ["N_CASES", "N_CONTROLS"]):
                if not None in [n_cases, n_controls]:
                    self.chunk_pl = self.chunk_pl.with_columns(
                    pl.lit(float(n_cases)).alias("N_CASES"),
                    pl.lit(float(n_controls)).alias("N_CONTROLS"),
                    pl.lit(float(n_cases) + float(n_controls)).alias("N"),
                    )
                else:
                    raise HarmonizationError("n_cases and n_controls columns are missing and were not provided")
                if self.mac is not None:
                    self.chunk_pl = self.chunk_pl.with_columns(
                    (2 * pl.col("N") * pl.min_horizontal(pl.col("EAF"), 1 - pl.col("EAF")))
                    .alias("MAC")
                    ).filter(pl.col("MAC") >= self.mac)
        else:
            raise HarmonizationError("Type of trait must be either binary or quant")


    
        if "SNPID" in self.chunk_pl.columns:
            self.chunk_pl = self.chunk_pl.drop("SNPID")
        
        
        #Here we always assumbe that the SNPs are in REF=A1 and ALT=A2
        swap = pl.col("A1") < pl.col("A2")

        #Start by creating new SNPID aligned
        self.chunk_pl = self.chunk_pl.with_columns([
        pl.when(swap).then(pl.col("A1")).otherwise(pl.col("A2")).alias("EA"),
        # set NEA to the larger allele
        pl.when(swap).then(pl.col("A2")).otherwise(pl.col("A1")).alias("NEA")]
        )

        self.chunk_pl = self.chunk_pl.with_columns(
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

        if self.pvar_file:
            self.align_alleles()
            self.chunk_pl.drop(["REF","ALT"])
        else:
            # flip the sign of BETA when swapping
            self.chunk_pl = self.chunk_pl.with_columns([
                pl.when(swap).then(-pl.col("BETA")).otherwise(pl.col("BETA")).alias("BETA"),
                # flip EAF to 1 - EAF when swapping
                pl.when(swap).then(1.0 - pl.col("EAF")).otherwise(pl.col("EAF")).alias("EAF"),
                # set EA to the smaller allele
            ])
        self.chunk_pl = self.chunk_pl.drop(["A1","A2"])

        self.chunk_pl = self.chunk_pl.with_columns(
                    pl.concat_str(
                        pl.lit("chr"),
                        pl.col("SNPID")).alias("SNPID"))
    
        if self.type_sumstat=="gwas":
            self.tiledb_types = {
                    "CHR": np.uint16,
                    "TRAIT": str,
                    "POS": np.uint32,
                    "SNPID": str,
                    "RSID": str,
                    "EAF": np.float32,
                    "BETA": np.float32,
                    "SE": np.float32,
                    "P": np.float64,
                }
            if "TRAIT" not in self.chunk_pl.columns:
                self.chunk_pl = self.chunk_pl.with_columns(
                    pl.lit(trait).alias("TRAIT")
                )
        else:
            self.tiledb_types = {
                "CHR": np.uint16,
                "CELL": str,
                "GENE": str,
                "POS": np.uint32,
                "SNPID": str,
                "RSID": str,
                "DIST": np.int64,
                "EAF": np.float32,
                "BETA": np.float32,
                "SE": np.float32,
                "P": np.float64,
            }
            if "CELL" not in self.chunk_pl.columns:
                self.chunk_pl = self.chunk_pl.with_columns(
                    pl.lit(cell).alias("CELL")
                )
            if "GENE" not in self.chunk_pl.columns:
                self.chunk_pl = self.chunk_pl.with_columns(
                    pl.lit(gene).alias("GENE")
                )
        if "RSID" not in self.chunk_pl.columns:
                self.chunk_pl = self.chunk_pl.with_columns(
                    pl.lit("None").alias("RSID")
                )
        if "LOG10P" in self.chunk_pl.columns:
            self.chunk_pl = self.chunk_pl.with_columns(
                    (10 ** (-pl.col("LOG10P"))).alias("P")
                )
        
        #Calculate p-value from z-score
        self.chunk_pl = self.chunk_pl.drop('P')
        self.chunk_pl = self.chunk_pl.with_columns(
            (pl.col("BETA") / pl.col("SE")).pow(2).map_batches(
            lambda x: pl.Series(stats.chi2.sf(x.to_numpy(), df=1)),
            return_dtype=pl.Float64
            ).alias('P')
            )
    
    def qc_sumstat(self, file_path:str):
        directory = self.uri + "_logs"
        filename = os.path.basename(file_path)   # "test.csv.gz"
        # Remove all extensions
        file_name = filename.split('.')[0]       # "test"
    
        if not os.path.isdir(directory):
            os.mkdir(directory)
        sumstat_preqc = self.chunk_pl.to_pandas()
        if self.type_sumstat == "gwas":
            if self.type_trait== "quant":
                sumstat_gl =gl.Sumstats(sumstat_preqc,
                    snpid="SNPID",
                    chrom="CHR",
                    pos="POS",
                    eaf="EAF",
                    beta="BETA",
                    se="SE",
                    p="P",
                    n="N",
                    ea = "EA",
                    nea = "NEA",
                    other = ["TRAIT","RSID"])
            else:
                sumstat_gl =gl.Sumstats(sumstat_preqc,
                    snpid="SNPID",
                    chrom="CHR",
                    pos="POS",
                    eaf="EAF",
                    beta="BETA",
                    se="SE",
                    p="P",
                    n="N",
                    ncase = "N_CASES",
                    ncontrol = "N_CONTROLS",
                    ea = "EA",
                    nea = "NEA",
                    other = ["TRAIT","RSID"])

        else:
            sumstat_gl =gl.Sumstats(sumstat_preqc,
                 snpid="SNPID",
                 chrom="CHR",
                 pos="POS",
                 eaf="EAF",
                 beta="BETA",
                 se="SE",
                 p="P",
                 n="N",
                 other = ["CELL","GENE","RSID","DIST","PHENO_VAR"])
        #sumstat_gl.fix_id()
        sumstat_gl.fix_chr(remove=True)
        sumstat_gl.fix_pos(remove=True)
        sumstat_gl.fix_allele(remove=True)
        sumstat_gl.check_sanity()
        sumstat_gl.check_data_consistency()
        #sumstat_gl.remove_dup(mode="m")
        #sumstat_gl.basic_check(n_cores = 4, remove=True, remove_dup=True)

        sumstat_gl.log.save(directory + "/" + file_name)
        self.chunk_pl = pl.from_pandas(sumstat_gl.data)
    
    def ingest_data(self, file_path):
        """Append harmonized data to TileDB."""
        pl.Config.set_tbl_cols(-1)
        if self.type_sumstat == "gwas":
            dedup_keys = ["CHR", "POS", "TRAIT"]
        else:
            dedup_keys = ["CHR", "POS", "CELL", "GENE"]

        # Option A (recommended): window count -> keep only rows whose group count == 1
        self.chunk_pl = (
            self.chunk_pl
                .with_columns(pl.count().over(dedup_keys).alias("_grp_count"))
                .filter(pl.col("_grp_count") == 1)
                .drop("_grp_count")
                )
        self.chunk_pl = self.chunk_pl.with_columns([
            pl.col("CHR").cast(pl.UInt16),
            pl.col("POS").cast(pl.UInt32)
        ])
        chunk_pl_ingest = self.chunk_pl.select(self.tiledb_types.keys())
        chunk_pl_ingest = chunk_pl_ingest.drop_nulls()
        try:
            tiledb.from_pandas(
                uri=self.uri,
                dataframe=chunk_pl_ingest.to_pandas(),
                index_dims=self.dimension_tiledb,
                column_types=self.tiledb_types,
                allows_duplicates = False,
                mode="append"
            )
            logger.info(f"Successfully appended chunk to TileDB for file {file_path}")
        except Exception as e:
            logger.error(f"Failed to append chunk to TileDB for file {file_path}: {e}")
            raise


    def create_metadata(self, file_path: str):
        """Create and store metadata as individual JSON files."""
        metadata = {
            "traits": [],
            "CELL": []
        }
        self.chunk_pl = self.chunk_pl.drop_nulls()
        if self.type_sumstat == "qtl":
            # Get unique cell types
            celltypes = self.chunk_pl["CELL"].unique().to_list()
            if not celltypes:
                raise HarmonizationError("No cell types found in the data")

            # Update CELL list
            metadata["CELL"] = celltypes

            # Process each cell type
            for cell in celltypes:
                # Filter by this cell type and compute ACAT per gene
                df_cell = self.chunk_pl.filter(pl.col("CELL") == cell)
            
                # Group by CHR and GENE, compute ACAT for each group
                chr_gene_agg = df_cell.group_by(["CHR", "GENE"]).agg([
                    pl.col("P").map_batches(
                        lambda s: pl.Series([acat_optimized(s)]),
                        return_dtype=pl.Float64
                    ).alias("ACAT_LIST"),
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
                    pheno_val = float(row["PHENO_VAR"])
                
                    if chrom not in metadata[cell]:
                        metadata[cell][chrom] = {}
                
                    gene_metadata = {
                        "ACAT": float(acat_val),
                        "N": n_val,
                        "PHENO_VAR": pheno_val
                    }
                
                    metadata[cell][chrom][gene] = gene_metadata
                    
        else:  # GWAS case
            if "TRAIT" not in self.chunk_pl.columns:
                raise HarmonizationError("TRAIT column is missing in the data")
        
            traits = self.chunk_pl["TRAIT"].unique().to_list()
            metadata["traits"] = traits
        
            for trait in traits:
                df_trait = self.chunk_pl.filter(pl.col("TRAIT") == trait)
                pheno_var = compute_pheno_variance(df_trait, self.type_trait)
                n_val = df_trait["N"].unique().to_list()[0]
            
                trait_metadata = {
                    "N": float(n_val),
                    "PHENO_VAR": float(pheno_var)
                }
            
                if self.type_trait == "binary":
                    n_cases = df_trait["N_CASES"].unique().to_list()[0]
                    n_controls = df_trait["N_CONTROLS"].unique().to_list()[0]
                    trait_metadata.update({
                        "N_CASES": float(n_cases),
                        "N_CONTROLS": float(n_controls)
                    })
            
                metadata[trait] = trait_metadata

        # Write individual metadata file instead of updating TileDB
        self._write_individual_metadata(metadata, file_path)
    
        logger.info(f"Individual metadata file created for {file_path}")

    def _write_individual_metadata(self, metadata, file_path):
        """Write metadata to individual JSON file."""
        metadata_dir = f"{self.uri}_metadata_parts"
        os.makedirs(metadata_dir, exist_ok=True)
    
        # Create a safe filename from the original file path
        file_stem = Path(file_path).stem
        metadata_file = os.path.join(metadata_dir, f"{file_stem}.json")
    
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
    
        logger.info(f"Metadata written to {metadata_file}")
    
    def export_metadata_to_csv(self, output_path: str = None):
        """Export metadata to CSV format."""
        if output_path is None:
            output_path = f"{self.uri}_metadata.csv"
    
        # Get the merged metadata from TileDB
        with tiledb.open(self.uri, "r") as array:
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
                                "PHENO_VAR": gene_metadata.get("PHENO_VAR", "")
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
                column_order = ["CHR", "CELL", "GENE", "ACAT", "N", "PHENO_VAR"]
                # Only include columns that exist in the DataFrame
                column_order = [col for col in column_order if col in df.columns]
                df = df[column_order]
            else:
                # GWAS format
                column_order = ["TRAIT", "N", "PHENO_VAR", "N_CASES", "N_CONTROLS"]
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

    def merge_metadata_files(self):
        """Merge all individual metadata files into final TileDB metadata."""
        metadata_dir = f"{self.uri}_metadata_parts"
    
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
        with tiledb.open(self.uri, "w") as array:
            array.meta["merged_metadata"] = json.dumps(merged_metadata)
    
        # Export to CSV
        self.export_metadata_to_csv()
    
        # Log summary
        cell_count = len(merged_metadata.get("CELL", []))
        trait_count = len(merged_metadata.get("traits", []))
    
        total_genes = 0
        for cell_type in merged_metadata.get("CELL", []):
            if cell_type in merged_metadata:
                for chrom_data in merged_metadata[cell_type].values():
                    total_genes += len(chrom_data)
    
        logger.info(f"Final merged metadata: {cell_count} cell types, {trait_count} traits, {total_genes} total genes")
    
