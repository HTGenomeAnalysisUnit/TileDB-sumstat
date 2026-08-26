import os 
import polars as pl
from tdbsumstat.utils.harmonize_ingest import Harmonize

def test_create_mapping(create_harmonized_obj):
    # 1. Execute the method
    create_harmonized_obj.create_mapping()
    # 2. Check the attribute on the object
    assert create_harmonized_obj.mapping_types == {
        "Chr": "CHR",
        "Gene": "GENE",
        "cell.type": "CELL",
        "pos": "POS",
        "a0": "A1",
        "a1": "A2",
        "p": "P",
        "N": "N",
        "beta": "BETA",
        "se": "SE",
    }

def test_create_tiledb(create_harmonized_obj, tmp_path):
    create_harmonized_obj.create_tiledb()
    # 2. Check the attribute on the object
    assert os.path.exists(f'{tmp_path}/test_tiledb_array')

def test_align_alleles(create_harmonized_obj, create_pvar):
    create_harmonized_obj.chunk_pl = pl.DataFrame({
        "CHR": 1,
        "POS": 1234,
        "SNPID": "1:1234:A:G",
        "BETA": -0.1, 
        "SE": 0.01,
        "EAF": 0.3
    })
    create_harmonized_obj.align_alleles()
    result_df = create_harmonized_obj.chunk_pl
    assert result_df.select("BETA").item() == 0.1
    assert result_df.select("EAF").item() == 0.7

