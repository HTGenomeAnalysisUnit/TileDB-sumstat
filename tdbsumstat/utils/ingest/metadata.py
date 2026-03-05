"""Metadata creation, export and merging for the ingestion pipeline."""
import json
import logging
import os
from pathlib import Path

import pandas as pd
import polars as pl
import tiledb

from tdbsumstat.utils import acat_optimized, compute_pheno_variance

logger = logging.getLogger(__name__)


class MetadataMixin:
    """Mixin that handles per-file and merged TileDB metadata."""

    # ------------------------------------------------------------------
    # Per-file metadata
    # ------------------------------------------------------------------

    def create_metadata(self, file_path: str) -> None:
        """Compute summary statistics and write a per-file JSON metadata fragment.

        The fragment is saved under ``{self.uri}_metadata_parts/``.

        Parameters
        ----------
        file_path:
            Original file path; used to derive the metadata filename.
        """
        from tdbsumstat.utils.ingest.errors import HarmonizationError

        metadata: dict = {"traits": [], "CELL": []}

        # Drop fully-null columns, then drop null rows
        self.chunk_pl = self.chunk_pl.select(
            [col for col in self.chunk_pl.columns if self.chunk_pl[col].null_count() < self.chunk_pl.height]
        )
        self.chunk_pl = self.chunk_pl.drop_nulls()

        if self.type_sumstat == "qtl":
            metadata = self._build_qtl_metadata(metadata)
        else:
            metadata = self._build_gwas_metadata(metadata)

        self._write_individual_metadata(metadata, file_path)
        logger.info("Individual metadata file created for %s", file_path)

    def _build_qtl_metadata(self, metadata: dict) -> dict:
        """Populate *metadata* with per-cell-type, per-gene statistics."""
        from tdbsumstat.utils.ingest.errors import HarmonizationError

        celltypes = self.chunk_pl["CELL"].unique().to_list()
        if not celltypes:
            raise HarmonizationError("No cell types found in the data")

        metadata["CELL"] = celltypes

        for cell in celltypes:
            df_cell = self.chunk_pl.filter(pl.col("CELL") == cell)

            chr_gene_agg = df_cell.group_by(["CHR", "GENE"]).agg([
                pl.col("P").map_batches(
                    lambda s: pl.Series([acat_optimized(s)]),
                    return_dtype=pl.Float64,
                ).alias("ACAT_LIST"),
                pl.col("P").min().alias("min_P"),
                pl.col("N").first().alias("N"),
                pl.col("PHENO_VAR").first().alias("PHENO_VAR"),
            ])
            chr_gene_agg = chr_gene_agg.with_columns(
                pl.col("ACAT_LIST").list.first().alias("ACAT")
            )

            if cell not in metadata:
                metadata[cell] = {}

            for row in chr_gene_agg.iter_rows(named=True):
                chrom = str(row["CHR"])
                gene = row["GENE"]
                metadata[cell].setdefault(chrom, {})
                metadata[cell][chrom][gene] = {
                    "ACAT": float(row["ACAT"]),
                    "N": float(row["N"]),
                    "PHENO_VAR": float(row["PHENO_VAR"]),
                    "MIN_P": float(row["min_P"]),
                }

        return metadata

    def _build_gwas_metadata(self, metadata: dict) -> dict:
        """Populate *metadata* with per-trait statistics."""
        from tdbsumstat.utils.ingest.errors import HarmonizationError

        if "TRAIT" not in self.chunk_pl.columns:
            raise HarmonizationError("TRAIT column is missing in the data")

        traits = self.chunk_pl["TRAIT"].unique().to_list()
        metadata["traits"] = traits

        for trait in traits:
            df_trait = self.chunk_pl.filter(pl.col("TRAIT") == trait)
            pheno_var = compute_pheno_variance(df_trait, self.type_trait)
            n_val = df_trait["N"].unique().to_list()[0]
            min_p = df_trait["P"].min().to_list()[0]

            trait_metadata: dict = {
                "N": float(n_val),
                "PHENO_VAR": float(pheno_var),
                "MIN_P": float(min_p),
            }

            if self.type_trait == "binary":
                trait_metadata["N_CASES"] = float(df_trait["N_CASES"].unique().to_list()[0])
                trait_metadata["N_CONTROLS"] = float(df_trait["N_CONTROLS"].unique().to_list()[0])

            metadata[trait] = trait_metadata

        return metadata

    def _write_individual_metadata(self, metadata: dict, file_path: str) -> None:
        """Persist a metadata dict as a JSON fragment under the metadata_parts dir."""
        metadata_dir = f"{self.uri}_metadata_parts"
        os.makedirs(metadata_dir, exist_ok=True)

        file_stem = Path(file_path).stem
        metadata_file = os.path.join(metadata_dir, f"{file_stem}.json")

        with open(metadata_file, "w") as fh:
            json.dump(metadata, fh, indent=2)

        logger.info("Metadata written to %s", metadata_file)

    # ------------------------------------------------------------------
    # Metadata export
    # ------------------------------------------------------------------

    def export_metadata_to_csv(self, output_path: str = None) -> pd.DataFrame:
        """Export the merged TileDB metadata to a CSV file.

        Parameters
        ----------
        output_path:
            Destination path.  Defaults to ``{self.uri}_metadata.csv``.

        Returns
        -------
        pd.DataFrame
            The exported metadata as a DataFrame (empty if nothing found).
        """
        if output_path is None:
            output_path = f"{self.uri}_metadata.csv"

        with tiledb.open(self.uri, "r") as array:
            merged_metadata_json = array.meta.get("merged_metadata", "{}")

        if not merged_metadata_json:
            logger.warning("No merged metadata found in TileDB array")
            return pd.DataFrame()

        merged_metadata = json.loads(merged_metadata_json)
        rows = []

        # QTL data
        for cell_type in merged_metadata.get("CELL", []):
            cell_data = merged_metadata.get(cell_type, {})
            for chrom, genes in cell_data.items():
                for gene_id, gm in genes.items():
                    rows.append({
                        "CHR": chrom, "CELL": cell_type, "GENE": gene_id,
                        "ACAT": gm.get("ACAT", ""),
                        "N": gm.get("N", ""),
                        "PHENO_VAR": gm.get("PHENO_VAR", ""),
                        "MIN_P": gm.get("MIN_P", ""),
                    })

        # GWAS data
        for trait in merged_metadata.get("traits", []):
            td = merged_metadata.get(trait, {})
            rows.append({
                "TRAIT": trait,
                "N": td.get("N", ""),
                "PHENO_VAR": td.get("PHENO_VAR", ""),
                "MIN_P": td.get("MIN_P", ""),
                "N_CASES": td.get("N_CASES", ""),
                "N_CONTROLS": td.get("N_CONTROLS", ""),
            })

        if not rows:
            logger.warning("No metadata found to export")
            return pd.DataFrame()

        df = pd.DataFrame(rows)

        if "CELL" in df.columns:
            col_order = [c for c in ["CHR", "CELL", "GENE", "ACAT", "MIN_P", "N", "PHENO_VAR"] if c in df.columns]
        else:
            col_order = [c for c in ["TRAIT", "N", "PHENO_VAR", "MIN_P", "N_CASES", "N_CONTROLS"] if c in df.columns]

        df = df[col_order]
        df.to_csv(output_path, index=False)
        logger.info("Metadata exported to %s", output_path)
        return df

    # ------------------------------------------------------------------
    # Merge all per-file fragments into final TileDB metadata
    # ------------------------------------------------------------------

    def merge_metadata_files(self) -> None:
        """Merge all per-file JSON fragments and store in TileDB.

        Reads every ``*.json`` from ``{self.uri}_metadata_parts/``,
        merges them into a single dict, stores it in ``TileDB.meta``
        under the key ``merged_metadata``, and also exports a CSV.
        """
        metadata_dir = f"{self.uri}_metadata_parts"

        if not os.path.exists(metadata_dir):
            logger.warning("No metadata directory found at %s", metadata_dir)
            return

        merged_metadata: dict = {"traits": [], "CELL": []}

        metadata_files = list(Path(metadata_dir).glob("*.json"))
        logger.info("Found %d metadata file(s) to merge", len(metadata_files))

        for metadata_file in metadata_files:
            try:
                with open(metadata_file, "r") as fh:
                    file_metadata = json.load(fh)
                self._merge_single_metadata(merged_metadata, file_metadata)
                logger.info("Processed %s", metadata_file.name)
            except Exception as exc:
                logger.error("Error processing metadata file %s: %s", metadata_file, exc)
                continue

        with tiledb.open(self.uri, "w") as array:
            array.meta["merged_metadata"] = json.dumps(merged_metadata)

        self.export_metadata_to_csv()

        cell_count = len(merged_metadata.get("CELL", []))
        trait_count = len(merged_metadata.get("traits", []))
        total_genes = sum(
            len(chrom_data)
            for cell_type in merged_metadata.get("CELL", [])
            for chrom_data in merged_metadata.get(cell_type, {}).values()
        )
        logger.info(
            "Final merged metadata: %d cell type(s), %d trait(s), %d total gene(s)",
            cell_count, trait_count, total_genes,
        )

    @staticmethod
    def _merge_single_metadata(merged: dict, file_meta: dict) -> None:
        """Merge one file_meta fragment into *merged* in-place."""
        # GWAS traits
        if file_meta.get("traits"):
            existing = set(merged.get("traits", []))
            merged["traits"] = list(existing.union(set(file_meta["traits"])))
            for trait in file_meta["traits"]:
                if trait in file_meta:
                    if trait in merged:
                        logger.warning("Trait %s already exists – overwriting", trait)
                    merged[trait] = file_meta[trait]

        # QTL cell types
        if file_meta.get("CELL"):
            existing = set(merged.get("CELL", []))
            merged["CELL"] = list(existing.union(set(file_meta["CELL"])))
            for cell in file_meta["CELL"]:
                if cell in file_meta:
                    merged.setdefault(cell, {})
                    for chrom, genes in file_meta[cell].items():
                        merged[cell].setdefault(chrom, {})
                        for gene, gene_data in genes.items():
                            if gene in merged[cell][chrom]:
                                logger.warning(
                                    "Gene %s in cell %s chr %s already exists – overwriting",
                                    gene, cell, chrom,
                                )
                            merged[cell][chrom][gene] = gene_data
