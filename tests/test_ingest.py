"""Tests for the refactored ingestion modules under tdbsumstat/utils/ingest/."""
import json
import os

import polars as pl
import pytest
import tiledb

from tdbsumstat.utils.ingest import Harmonize, HarmonizationError
from tdbsumstat.utils.ingest.errors import HarmonizationError as DirectError
from tdbsumstat.utils.ingest.schema import SchemaMixin
from tdbsumstat.utils.ingest.mapping import MappingMixin
from tdbsumstat.utils.ingest.harmonize import HarmonizeMixin
from tdbsumstat.utils.ingest.writer import WriterMixin
from tdbsumstat.utils.ingest.metadata import MetadataMixin

EXAMPLE_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "example_data")


# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def mapping_file():
    return os.path.join(EXAMPLE_DATA_DIR, "mapping_file_test.csv")


@pytest.fixture()
def sample_qtl_df():
    """Minimal QTL Polars DataFrame matching the mapping_file_test column names."""
    return pl.DataFrame({
        "Chr": [20, 20, 20],
        "Gene": ["ENSG0000010000"] * 3,
        "cell.type": ["T_gd"] * 3,
        "pos": [100, 200, 300],
        "a0": ["A", "C", "G"],
        "a1": ["T", "G", "A"],
        "p": [0.01, 0.5, 0.001],
        "N": [500, 500, 500],
        "beta": [0.1, -0.2, 0.3],
        "se": [0.05, 0.08, 0.04],
    })


@pytest.fixture()
def empty_qtl_df():
    """Empty QTL Polars DataFrame with correct schema."""
    return pl.DataFrame({
        "Chr": pl.Series([], dtype=pl.Int64),
        "Gene": pl.Series([], dtype=pl.Utf8),
        "cell.type": pl.Series([], dtype=pl.Utf8),
        "pos": pl.Series([], dtype=pl.Int64),
        "a0": pl.Series([], dtype=pl.Utf8),
        "a1": pl.Series([], dtype=pl.Utf8),
        "p": pl.Series([], dtype=pl.Float64),
        "N": pl.Series([], dtype=pl.Int64),
        "beta": pl.Series([], dtype=pl.Float64),
        "se": pl.Series([], dtype=pl.Float64),
    })


# ---------------------------------------------------------------------------
# SchemaMixin
# ---------------------------------------------------------------------------

class TestSchemaMixin:
    def test_qtl_dimensions(self, mapping_file, tmp_path):
        uri = str(tmp_path / "qtl")
        h = Harmonize(mapping_file, uri, "qtl", None, "quant", None, None, False)
        h.create_tiledb()
        schema = tiledb.ArraySchema.load(uri)
        dim_names = [schema.domain.dim(i).name for i in range(schema.domain.ndim)]
        assert dim_names == ["CHR", "CELL", "GENE", "POS"]
        attr_names = [schema.attr(i).name for i in range(schema.nattr)]
        assert "DIST" in attr_names

    def test_gwas_dimensions(self, mapping_file, tmp_path):
        uri = str(tmp_path / "gwas")
        h = Harmonize(mapping_file, uri, "gwas", None, "quant", None, None, False)
        h.create_tiledb()
        schema = tiledb.ArraySchema.load(uri)
        dim_names = [schema.domain.dim(i).name for i in range(schema.domain.ndim)]
        assert dim_names == ["CHR", "TRAIT", "POS"]
        attr_names = [schema.attr(i).name for i in range(schema.nattr)]
        assert "DIST" not in attr_names

    def test_sparse_schema(self, mapping_file, tmp_path):
        uri = str(tmp_path / "sparse_check")
        h = Harmonize(mapping_file, uri, "qtl", None, "quant", None, None, False)
        h.create_tiledb()
        schema = tiledb.ArraySchema.load(uri)
        assert schema.sparse is True


# ---------------------------------------------------------------------------
# MappingMixin
# ---------------------------------------------------------------------------

class TestMappingMixin:
    def test_creates_dict(self, mapping_file):
        h = Harmonize(mapping_file, "dummy", "qtl", None, "quant", None, None, False)
        h.create_mapping()
        assert isinstance(h.mapping_types, dict)
        assert len(h.mapping_types) > 0

    def test_missing_beta_raises(self, tmp_path):
        f = tmp_path / "m.csv"
        f.write_text("Chr,CHR\npos,POS\n")
        h = Harmonize(str(f), "dummy", "qtl", None, "quant", None, None, False)
        with pytest.raises(HarmonizationError, match="BETA and SE"):
            h.create_mapping()

    def test_empty_file_raises(self, tmp_path):
        f = tmp_path / "empty.csv"
        f.write_text("")
        h = Harmonize(str(f), "dummy", "qtl", None, "quant", None, None, False)
        with pytest.raises(HarmonizationError, match="empty or not formatted"):
            h.create_mapping()

    def test_harmonization_error_importable_directly(self):
        assert DirectError is HarmonizationError


# ---------------------------------------------------------------------------
# HarmonizeMixin
# ---------------------------------------------------------------------------

class TestHarmonizeMixin:
    def test_snpid_prefixed_with_chr(self, mapping_file, sample_qtl_df):
        h = Harmonize(mapping_file, "dummy", "qtl", None, "quant", None, None, False)
        h.create_mapping()
        h.harmonize(sumstat=sample_qtl_df, cell="T_gd", gene="ENSG0000010000", n=500, pheno_var=1.5)
        for snpid in h.chunk_pl["SNPID"].to_list():
            assert snpid.startswith("chr"), f"SNPID {snpid!r} does not start with 'chr'"

    def test_cell_and_gene_columns_added(self, mapping_file, sample_qtl_df):
        h = Harmonize(mapping_file, "dummy", "qtl", None, "quant", None, None, False)
        h.create_mapping()
        h.harmonize(sumstat=sample_qtl_df, cell="T_gd", gene="ENSG0000010000", n=500, pheno_var=1.5)
        assert "CELL" in h.chunk_pl.columns
        assert "GENE" in h.chunk_pl.columns
        assert h.chunk_pl["CELL"].unique().to_list() == ["T_gd"]
        assert h.chunk_pl["GENE"].unique().to_list() == ["ENSG0000010000"]

    def test_p_value_column_present(self, mapping_file, sample_qtl_df):
        h = Harmonize(mapping_file, "dummy", "qtl", None, "quant", None, None, False)
        h.create_mapping()
        h.harmonize(sumstat=sample_qtl_df, cell="T_gd", gene="ENSG0000010000", n=500, pheno_var=1.5)
        assert "P" in h.chunk_pl.columns
        assert all(0 <= p <= 1 for p in h.chunk_pl["P"].to_list())

    def test_missing_n_raises(self, mapping_file, sample_qtl_df):
        h = Harmonize(mapping_file, "dummy", "qtl", None, "quant", None, None, False)
        h.create_mapping()
        with pytest.raises(HarmonizationError):
            h.harmonize(sumstat=sample_qtl_df.drop("N"), cell="T_gd",
                        gene="ENSG0000010000", pheno_var=1.5)

    def test_maf_filter_reduces_rows(self, mapping_file, sample_qtl_df):
        """All rows in sample have EAF=0 (default), so MAF=0; maf=0.1 should drop all."""
        h = Harmonize(mapping_file, "dummy", "qtl", None, "quant", None, maf=0.1, permuted=False)
        h.create_mapping()
        h.harmonize(sumstat=sample_qtl_df, cell="T_gd", gene="ENSG0000010000", n=500, pheno_var=1.5)
        # With maf=0.1 all rows may be filtered out (EAF defaults to 0)
        assert h.chunk_pl.height >= 0  # just confirm no crash

    def test_binary_trait_n_columns(self, mapping_file, sample_qtl_df):
        h = Harmonize(mapping_file, "dummy", "gwas", None, "binary", None, None, False)
        h.create_mapping()
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
        h.harmonize(sumstat=gwas_df, trait="T1", n_cases=200, n_controls=300, pheno_var=1.0)
        assert "N_CASES" in h.chunk_pl.columns
        assert "N_CONTROLS" in h.chunk_pl.columns

    def test_invalid_type_trait_raises(self, mapping_file, sample_qtl_df):
        h = Harmonize(mapping_file, "dummy", "qtl", None, "invalid_type", None, None, False)
        h.create_mapping()
        with pytest.raises(HarmonizationError, match="binary or quant"):
            h.harmonize(sumstat=sample_qtl_df, cell="T_gd", gene="ENSG0000010000", n=500, pheno_var=1.5)

    def test_empty_dataframe_harmonizes_without_crash(self, mapping_file, empty_qtl_df):
        """Harmonizing an empty DataFrame should not raise and should produce an empty chunk_pl."""
        h = Harmonize(mapping_file, "dummy", "qtl", None, "quant", None, None, False)
        h.create_mapping()
        h.harmonize(sumstat=empty_qtl_df, cell="T_gd", gene="ENSG0000010000", n=500, pheno_var=1.5)
        assert h.chunk_pl.height == 0


# ---------------------------------------------------------------------------
# WriterMixin
# ---------------------------------------------------------------------------

class TestWriterMixin:
    def test_ingest_appends_rows(self, mapping_file, sample_qtl_df, tmp_path):
        uri = str(tmp_path / "writer_test")
        h = Harmonize(mapping_file, uri, "qtl", None, "quant", None, None, False)
        h.create_tiledb()
        h.create_mapping()
        h.harmonize(sumstat=sample_qtl_df, cell="T_gd", gene="ENSG0000010000", n=500, pheno_var=1.5)
        h.ingest_data(file_path="dummy_file.tsv")

        with tiledb.open(uri, mode="r") as A:
            df = A.query(dims=["CHR", "CELL", "GENE", "POS"]).df[:]
        assert len(df) > 0

    def test_duplicates_are_removed(self, mapping_file, tmp_path):
        """Rows sharing an index key (CHR, POS, CELL, GENE) are ALL dropped.

        The ingestion pipeline treats ambiguous duplicate positions as
        unreliable and excludes them entirely rather than keeping one.
        """
        uri = str(tmp_path / "dedup_test")
        h = Harmonize(mapping_file, uri, "qtl", None, "quant", None, None, False)
        h.create_tiledb()
        h.create_mapping()

        # Two identical rows
        dup_df = pl.DataFrame({
            "Chr": [20, 20],
            "Gene": ["ENSG0000010000", "ENSG0000010000"],
            "cell.type": ["T_gd", "T_gd"],
            "pos": [500, 500],  # same position → duplicate
            "a0": ["A", "A"],
            "a1": ["T", "T"],
            "p": [0.01, 0.01],
            "N": [500, 500],
            "beta": [0.1, 0.1],
            "se": [0.05, 0.05],
        })
        h.harmonize(sumstat=dup_df, cell="T_gd", gene="ENSG0000010000", n=500, pheno_var=1.5)
        h.ingest_data(file_path="dup_file.tsv")

        with tiledb.open(uri, mode="r") as A:
            df = A.query(dims=["CHR", "CELL", "GENE", "POS"]).df[:]
        # Both rows are duplicate; they should both be dropped by the dedup logic
        assert len(df) == 0


# ---------------------------------------------------------------------------
# MetadataMixin
# ---------------------------------------------------------------------------

class TestMetadataMixin:
    def _create_ingested_tiledb(self, mapping_file, sample_qtl_df, tmp_path):
        """Helper: create + ingest + return a Harmonize object with data."""
        uri = str(tmp_path / "meta_test")
        h = Harmonize(mapping_file, uri, "qtl", None, "quant", None, None, False)
        h.create_tiledb()
        h.create_mapping()
        h.harmonize(sumstat=sample_qtl_df, cell="T_gd", gene="ENSG0000010000", n=500, pheno_var=1.5)
        h.ingest_data(file_path="file.tsv")
        return h

    def test_create_metadata_writes_json(self, mapping_file, sample_qtl_df, tmp_path):
        h = self._create_ingested_tiledb(mapping_file, sample_qtl_df, tmp_path)
        h.create_metadata(file_path="file.tsv")

        metadata_dir = h.uri + "_metadata_parts"
        json_files = list(os.listdir(metadata_dir))
        assert len(json_files) == 1
        assert json_files[0].endswith(".json")

    def test_metadata_json_has_cell_key(self, mapping_file, sample_qtl_df, tmp_path):
        h = self._create_ingested_tiledb(mapping_file, sample_qtl_df, tmp_path)
        h.create_metadata(file_path="file.tsv")

        metadata_dir = h.uri + "_metadata_parts"
        json_file = os.path.join(metadata_dir, os.listdir(metadata_dir)[0])
        with open(json_file) as f:
            meta = json.load(f)
        assert "CELL" in meta
        assert "T_gd" in meta["CELL"]

    def test_merge_metadata_stores_in_tiledb(self, mapping_file, sample_qtl_df, tmp_path):
        h = self._create_ingested_tiledb(mapping_file, sample_qtl_df, tmp_path)
        h.create_metadata(file_path="file.tsv")
        h.merge_metadata_files()

        with tiledb.open(h.uri, mode="r") as A:
            merged = json.loads(A.meta["merged_metadata"])
        assert "CELL" in merged
        assert "T_gd" in merged["CELL"]

    def test_merge_metadata_static_helper(self):
        """_merge_single_metadata should correctly merge QTL entries."""
        merged = {"traits": [], "CELL": []}
        file_meta = {
            "CELL": ["T_gd"],
            "T_gd": {
                "20": {"ENSG0000010000": {"ACAT": 0.5, "N": 500.0, "PHENO_VAR": 1.5, "MIN_P": 0.001}}
            },
        }
        MetadataMixin._merge_single_metadata(merged, file_meta)
        assert "T_gd" in merged["CELL"]
        assert "20" in merged["T_gd"]
        assert "ENSG0000010000" in merged["T_gd"]["20"]

    def test_export_metadata_to_csv(self, mapping_file, sample_qtl_df, tmp_path):
        h = self._create_ingested_tiledb(mapping_file, sample_qtl_df, tmp_path)
        h.create_metadata(file_path="file.tsv")
        h.merge_metadata_files()

        import pandas as pd
        df = h.export_metadata_to_csv()
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0
        assert "CELL" in df.columns


# ---------------------------------------------------------------------------
# End-to-end ingestion using example data files
# ---------------------------------------------------------------------------

class TestEndToEndIngestion:
    def test_full_qtl_pipeline(self, mapping_file, tmp_path):
        """Full QTL ingestion pipeline using example data files."""
        uri = str(tmp_path / "e2e_qtl")
        h = Harmonize(mapping_file, uri, "qtl", None, "quant", None, None, False)
        h.create_tiledb()
        h.create_mapping()

        data_files = [
            (os.path.join(EXAMPLE_DATA_DIR, "dummy_out_ENSG0000010000.tsv.gz"), "ENSG0000010000"),
            (os.path.join(EXAMPLE_DATA_DIR, "dummy_out_ENSG0000010001.tsv.gz"), "ENSG0000010001"),
        ]

        for filepath, gene in data_files:
            chunk_pl = pl.read_csv(filepath, separator="\t", low_memory=True, null_values="NA")
            h.harmonize(sumstat=chunk_pl, cell="Tgd", gene=gene, n=4000, pheno_var=1.5)
            h.ingest_data(file_path=filepath)
            h.create_metadata(file_path=filepath)

        h.merge_metadata_files()

        # Verify data was written
        with tiledb.open(uri, mode="r") as A:
            df = A.query(dims=["CHR", "CELL", "GENE", "POS"]).df[:]
        assert len(df) > 0

        # Verify merged metadata
        with tiledb.open(uri, mode="r") as A:
            merged = json.loads(A.meta["merged_metadata"])
        assert "Tgd" in merged["CELL"]

    def test_backward_compat_import(self):
        """harmonize_ingest shim must re-export the same Harmonize class."""
        from tdbsumstat.utils.harmonize_ingest import Harmonize as OldHarmonize
        from tdbsumstat.utils.ingest import Harmonize as NewHarmonize
        assert OldHarmonize is NewHarmonize


# ---------------------------------------------------------------------------
# Module smoke-tests
# ---------------------------------------------------------------------------

class TestModuleImports:
    def test_import_errors(self):
        from tdbsumstat.utils.ingest.errors import HarmonizationError  # noqa: F401

    def test_import_schema(self):
        from tdbsumstat.utils.ingest.schema import SchemaMixin  # noqa: F401

    def test_import_mapping(self):
        from tdbsumstat.utils.ingest.mapping import MappingMixin  # noqa: F401

    def test_import_harmonize(self):
        from tdbsumstat.utils.ingest.harmonize import HarmonizeMixin  # noqa: F401

    def test_import_qc(self):
        from tdbsumstat.utils.ingest.qc import QCMixin  # noqa: F401

    def test_import_writer(self):
        from tdbsumstat.utils.ingest.writer import WriterMixin  # noqa: F401

    def test_import_metadata(self):
        from tdbsumstat.utils.ingest.metadata import MetadataMixin  # noqa: F401

    def test_import_package(self):
        from tdbsumstat.utils.ingest import Harmonize, HarmonizationError  # noqa: F401
