"""Pytest fixtures shared across the test suite."""
import os
import tempfile

import numpy as np
import pandas as pd
import polars as pl
import pytest
import tiledb

EXAMPLE_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "example_data")


@pytest.fixture(scope="session")
def example_data_dir():
    """Return the path to the example_data directory."""
    return os.path.abspath(EXAMPLE_DATA_DIR)


@pytest.fixture(scope="session")
def mapping_file_sc(example_data_dir):
    """Return the path to the QTL mapping file."""
    return os.path.join(example_data_dir, "mapping_file_test.csv")


@pytest.fixture(scope="session")
def qtl_tiledb(tmp_path_factory, example_data_dir, mapping_file_sc):
    """Create a QTL TileDB from example data and return its path.

    This fixture runs once per test session and reuses the same TileDB.
    """
    from tdbsumstat.utils.harmonize_ingest import Harmonize

    uri = str(tmp_path_factory.mktemp("tiledb") / "test_qtl_tiledb")

    h = Harmonize(
        mapping_file=mapping_file_sc,
        uri=uri,
        type_sumstat="qtl",
        pvar_file=None,
        type_trait="quant",
        mac=None,
        maf=None,
        permuted=False,
    )
    h.create_tiledb()
    h.create_mapping()

    # Ingest both example files
    data_files = [
        os.path.join(example_data_dir, "dummy_out_ENSG0000010000.tsv.gz"),
        os.path.join(example_data_dir, "dummy_out_ENSG0000010001.tsv.gz"),
    ]

    for filepath in data_files:
        chunk_pl = pl.read_csv(filepath, separator="\t", low_memory=True, null_values="NA")
        # Extract gene from filename: dummy_out_ENSG0000010000.tsv.gz -> ENSG0000010000
        gene = os.path.basename(filepath).replace("dummy_out_", "").split(".")[0]
        h.harmonize(sumstat=chunk_pl, cell="T_gd", gene=gene, n=4000, pheno_var=1.5)
        h.ingest_data(file_path=filepath)
        h.create_metadata(file_path=filepath)

    h.merge_metadata_files()
    return uri


@pytest.fixture(scope="session")
def qtl_metadata_df(qtl_tiledb):
    """Return the metadata DataFrame for the QTL TileDB."""
    import json

    with tiledb.open(qtl_tiledb, mode="r") as A:
        metadata = json.loads(A.meta["merged_metadata"])

    rows = []
    for cell in metadata["CELL"]:
        for chrom, genes in metadata[cell].items():
            for gene, stats in genes.items():
                rows.append({"CELL": cell, "CHR": chrom, "GENE": gene, **stats})
    df_meta = pl.DataFrame(rows)
    df_meta = df_meta.with_columns(pl.col("CHR").cast(pl.UInt16))
    return df_meta


@pytest.fixture()
def tmp_out(tmp_path):
    """Return a temporary output prefix for test files."""
    return str(tmp_path / "out")
