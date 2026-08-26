import pytest
import pandas as pd
from tdbsumstat.utils.harmonize_ingest import Harmonize

@pytest.fixture
def create_df_sumstat():
    ensg00 = pd.read_csv('tests/data/dummy_out_ENSG0000010001.tsv.gz', sep = '\t')
    return ensg00

@pytest.fixture
def create_harmonized_obj(tmp_path):
    obj = Harmonize(
        mapping_file= 'tests/data/mapping_file_test.csv',
        uri = str(tmp_path / "test_tiledb_array"),
        type_sumstat = "qtl",
        pvar_file = str(tmp_path / "temp_pvar.pvar"),
        type_trait = 'quant',
        permuted = False,
        mac = 10,
        maf = 0.001,
            )
    return obj

@pytest.fixture
def create_pvar(tmp_path):
    # Added list brackets to ensure pandas creates rows correctly
    pvar = pd.DataFrame({
        "CHROM": ["1"], 
        "POS": ["1234"], 
        "SNPID": ["1:1234:A:G"], 
        "REF": ["A"],
        "ALT": ["G"]
    })
    
    # CRITICAL: index=False prevents pandas from writing a row-number column
    pvar_path = f"{tmp_path}/temp_pvar.pvar"
    pvar.to_csv(pvar_path, sep='\t', index=False)
    
    return pvar_path