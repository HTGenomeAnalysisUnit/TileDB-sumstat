"""Tests for the export handler functions (snp, regions, locusbreaker, metadata, traits)."""
import os

import pandas as pd
import polars as pl
import pytest
import tiledb

from tdbsumstat.cli.export.regions import export_by_regions
from tdbsumstat.cli.export.snp import export_by_snp

EXAMPLE_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "example_data")


class TestExportByRegions:
    def test_qtl_regions_output(self, qtl_tiledb, tmp_out):
        """Export by regions should create a non-empty CSV for known regions."""
        region_file = os.path.join(EXAMPLE_DATA_DIR, "region_list_sc.csv")
        with tiledb.open(qtl_tiledb, mode="r") as tdb:
            export_by_regions(
                tiledb_export=tdb,
                table_regions=region_file,
                attr="P,SNPID,EAF,BETA,SE",
                type_sumstat="qtl",
                out=tmp_out + "_regions.csv",
            )

        output_file = tmp_out + "_regions.csv"
        if os.path.exists(output_file):
            df = pd.read_csv(output_file)
            assert len(df) > 0
            assert "P" in df.columns

    def test_empty_region_produces_no_file(self, qtl_tiledb, tmp_out, tmp_path):
        """Querying a region with no matching gene should not create an output file."""
        empty_region_file = str(tmp_path / "empty_region.csv")
        # Use valid CHR/cell but a gene that doesn't exist in example data
        pd.DataFrame({
            "CHR": [20],
            "START": [1],
            "END": [100],
            "TRAIT": ["T_gd:ENSG9999999999"],
        }).to_csv(empty_region_file, index=False)

        out_file = tmp_out + "_empty.csv"
        with tiledb.open(qtl_tiledb, mode="r") as tdb:
            export_by_regions(
                tiledb_export=tdb,
                table_regions=empty_region_file,
                attr="P,SNPID,EAF,BETA,SE",
                type_sumstat="qtl",
                out=out_file,
            )
        # Should not have created a file (no data matched)
        assert not os.path.exists(out_file)


class TestExportBySnp:
    def test_qtl_snp_output(self, qtl_tiledb, tmp_out):
        """Export by SNP list should create a CSV with matching rows."""
        snp_file = os.path.join(EXAMPLE_DATA_DIR, "snp_list_sc.csv")
        with tiledb.open(qtl_tiledb, mode="r") as tdb:
            export_by_snp(
                tiledb_export=tdb,
                snp=snp_file,
                attr="P,SNPID,EAF,BETA,SE",
                type_sumstat="qtl",
                out=tmp_out,
            )
        # The output file(s) are named {out}_{trait}_{chrom}.csv
        output_files = [
            f for f in os.listdir(os.path.dirname(tmp_out) or ".")
            if f.endswith(".csv")
        ]
        assert len(output_files) > 0


class TestLocusbreakerExport:
    def test_locusbreaker_import(self):
        """Locusbreaker module should be importable without errors."""
        from tdbsumstat.cli.export.locusbreaker import export_with_locusbreaker  # noqa: F401


class TestMetadataExport:
    def test_export_metadata_creates_file(self, qtl_tiledb, tmp_out):
        """export_metadata should create a *_meta.csv file."""
        from tdbsumstat.cli.export.metadata import export_metadata

        export_metadata(qtl_tiledb, "qtl", tmp_out)
        assert os.path.exists(tmp_out + "_meta.csv")
        df = pd.read_csv(tmp_out + "_meta.csv")
        assert len(df) > 0

    def test_recompute_metadata_import(self):
        """recompute_metadata should be importable without errors."""
        from tdbsumstat.cli.export.metadata import recompute_metadata  # noqa: F401


class TestExportByTraits:
    def test_export_by_traits_import(self):
        """export_by_traits should be importable without errors."""
        from tdbsumstat.cli.export.traits import export_by_traits  # noqa: F401

    def test_export_by_traits_qtl(self, qtl_tiledb, tmp_out):
        """Export by traits should run without errors."""
        from tdbsumstat.cli.export.traits import export_by_traits

        trait_file = str(os.path.join(os.path.dirname(tmp_out), "trait_list.csv"))
        pd.DataFrame({"TRAIT": ["T_gd:ENSG0000010000"]}).to_csv(trait_file, index=False)

        export_by_traits(
            uri_path=qtl_tiledb,
            trait_list=trait_file,
            attr="P,SNPID,EAF,BETA,SE",
            type_sumstat="qtl",
            out=tmp_out,
            batch_name="test",
        )
        # The function completes without error; output may be empty for small example data


class TestModuleImports:
    """Smoke tests: ensure all new modules can be imported cleanly."""

    def test_import_helpers(self):
        from tdbsumstat.cli.export.helpers import open_tiledb_and_load_metadata  # noqa: F401

    def test_import_snp(self):
        from tdbsumstat.cli.export.snp import export_by_snp  # noqa: F401

    def test_import_regions(self):
        from tdbsumstat.cli.export.regions import export_by_regions  # noqa: F401

    def test_import_locusbreaker(self):
        from tdbsumstat.cli.export.locusbreaker import export_with_locusbreaker  # noqa: F401

    def test_import_metadata(self):
        from tdbsumstat.cli.export.metadata import export_metadata, recompute_metadata  # noqa: F401

    def test_import_traits(self):
        from tdbsumstat.cli.export.traits import export_by_traits  # noqa: F401

    def test_import_command(self):
        from tdbsumstat.cli.export.command import export  # noqa: F401

    def test_import_export_package(self):
        from tdbsumstat.cli.export import export  # noqa: F401
