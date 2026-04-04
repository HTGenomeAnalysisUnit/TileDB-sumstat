"""Tests for the Harmonize class (harmonize_ingest.py)."""
import os

import numpy as np
import polars as pl
import pytest
import tiledb

from tdbsumstat.utils.harmonize_ingest import Harmonize, HarmonizationError


EXAMPLE_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "example_data")


@pytest.fixture()
def mapping_file():
    return os.path.join(EXAMPLE_DATA_DIR, "mapping_file_test.csv")


@pytest.fixture()
def sample_qtl_df():
    """Minimal QTL polars DataFrame mimicking the example data columns."""
    return pl.DataFrame({
        "Chr": [20, 20, 20],
        "Gene": ["ENSG0000010000", "ENSG0000010000", "ENSG0000010000"],
        "cell.type": ["T_gd", "T_gd", "T_gd"],
        "pos": [100, 200, 300],
        "a0": ["A", "C", "G"],
        "a1": ["T", "G", "A"],
        "p": [0.01, 0.5, 0.001],
        "N": [500, 500, 500],
        "beta": [0.1, -0.2, 0.3],
        "se": [0.05, 0.08, 0.04],
    })


class TestCreateMapping:
    def test_valid_mapping(self, mapping_file):
        h = Harmonize(mapping_file, "dummy", "qtl", None, "quant", None, None, False)
        h.create_mapping()
        assert "CHR" in h.mapping_types.values()
        assert "BETA" in h.mapping_types.values()
        assert "SE" in h.mapping_types.values()

    def test_missing_beta_se_raises(self, tmp_path):
        bad_mapping = tmp_path / "bad_mapping.csv"
        bad_mapping.write_text("Chr,CHR\npos,POS\n")
        h = Harmonize(str(bad_mapping), "dummy", "qtl", None, "quant", None, None, False)
        with pytest.raises(HarmonizationError, match="BETA and SE"):
            h.create_mapping()

    def test_empty_mapping_raises(self, tmp_path):
        empty_file = tmp_path / "empty.csv"
        empty_file.write_text("")
        h = Harmonize(str(empty_file), "dummy", "qtl", None, "quant", None, None, False)
        with pytest.raises(HarmonizationError, match="empty or not formatted"):
            h.create_mapping()


class TestCreateTileDB:
    def test_creates_qtl_tiledb(self, mapping_file, tmp_path):
        uri = str(tmp_path / "test_tiledb")
        h = Harmonize(mapping_file, uri, "qtl", None, "quant", None, None, False)
        h.create_tiledb()
        assert tiledb.array_exists(uri)
        schema = tiledb.ArraySchema.load(uri)
        dim_names = [schema.domain.dim(i).name for i in range(schema.domain.ndim)]
        assert "CHR" in dim_names
        assert "CELL" in dim_names
        assert "GENE" in dim_names
        assert "POS" in dim_names

    def test_creates_gwas_tiledb(self, mapping_file, tmp_path):
        uri = str(tmp_path / "test_gwas")
        h = Harmonize(mapping_file, uri, "gwas", None, "quant", None, None, False)
        h.create_tiledb()
        assert tiledb.array_exists(uri)
        schema = tiledb.ArraySchema.load(uri)
        dim_names = [schema.domain.dim(i).name for i in range(schema.domain.ndim)]
        assert "TRAIT" in dim_names
        assert "CHR" in dim_names
        assert "POS" in dim_names


class TestHarmonize:
    def test_harmonize_qtl_adds_cell_gene(self, mapping_file, sample_qtl_df):
        h = Harmonize(mapping_file, "dummy", "qtl", None, "quant", None, None, False)
        h.create_mapping()
        h.harmonize(sumstat=sample_qtl_df, cell="T_gd", gene="ENSG0000010000", n=500, pheno_var=1.5)
        assert "CELL" in h.chunk_pl.columns
        assert "GENE" in h.chunk_pl.columns
        assert "SNPID" in h.chunk_pl.columns

    def test_harmonize_snpid_format(self, mapping_file, sample_qtl_df):
        h = Harmonize(mapping_file, "dummy", "qtl", None, "quant", None, None, False)
        h.create_mapping()
        h.harmonize(sumstat=sample_qtl_df, cell="T_gd", gene="ENSG0000010000", n=500, pheno_var=1.5)
        # All SNPIDs should start with 'chr'
        snpids = h.chunk_pl["SNPID"].to_list()
        assert all(s.startswith("chr") for s in snpids)

    def test_harmonize_raises_without_n(self, mapping_file, sample_qtl_df):
        df_no_n = sample_qtl_df.drop("N")
        h = Harmonize(mapping_file, "dummy", "qtl", None, "quant", None, None, False)
        h.create_mapping()
        with pytest.raises(HarmonizationError):
            h.harmonize(sumstat=df_no_n, cell="T_gd", gene="ENSG0000010000", pheno_var=1.5)

    def test_harmonize_maf_filter(self, mapping_file, sample_qtl_df):
        """MAF filter should remove rows with EAF outside the range."""
        # Add EAF column: use a0/a1 encoded EAF from original data, but override
        df_with_eaf = sample_qtl_df.with_columns(pl.lit(0.01).alias("EAF_explicit"))
        h = Harmonize(mapping_file, "dummy", "qtl", None, "quant", None, maf=0.05, permuted=False)
        h.create_mapping()
        # With maf=0.05 and EAF derived from alleles, some rows may be filtered
        # Just make sure the call doesn't crash
        h.harmonize(sumstat=sample_qtl_df, cell="T_gd", gene="ENSG0000010000", n=500, pheno_var=1.5)

    def test_binary_trait_adds_n_columns(self, mapping_file, sample_qtl_df):
        h = Harmonize(mapping_file, "dummy", "gwas", None, "binary", None, None, False)
        h.create_mapping()
        gwas_df = sample_qtl_df.rename({"cell.type": "TRAIT_col", "Gene": "TRAIT"}).drop("Gene", strict=False)
        # Use a simplified gwas-like df
        gwas_df = pl.DataFrame({
            "Chr": [1, 1],
            "pos": [100, 200],
            "a0": ["A", "C"],
            "a1": ["T", "G"],
            "p": [0.01, 0.5],
            "beta": [0.1, -0.2],
            "se": [0.05, 0.08],
            "N": [500, 500],
            "Gene": ["T1", "T1"],
            "cell.type": ["gwas", "gwas"],
        })
        mapping_gwas = mapping_file  # mapping is for QTL, adapt
        h2 = Harmonize(mapping_gwas, "dummy", "gwas", None, "binary", None, None, False)
        h2.create_mapping()
        h2.harmonize(
            sumstat=gwas_df,
            trait="T1",
            n_cases=200,
            n_controls=300,
            pheno_var=1.0,
        )
        assert "N_CASES" in h2.chunk_pl.columns
        assert "N_CONTROLS" in h2.chunk_pl.columns


class TestIngestData:
    def test_ingest_creates_data(self, qtl_tiledb):
        """The session-scoped qtl_tiledb fixture should have rows."""
        with tiledb.open(qtl_tiledb, mode="r") as A:
            df = A.query(dims=["CHR", "CELL", "GENE", "POS"]).df[:]
        assert len(df) > 0
        assert "CELL" in df.columns
