# TileDB for summary statistics analysis

<img width="345" alt="Screenshot 2024-11-17 at 20 47 08" src="https://github.com/user-attachments/assets/cd7ffb0e-1d58-41a9-9b86-d29420849111">

### This is a program for the and analysis of single cell QTL and GWAS data using Nextflow and TileDB

#### Description

Welcome to TileDB-sumstat. This is a custom TileDB for singel cell QTL and GWAS data.
Please have a look [here](https://github.com/HTGenomeAnalysisUnit/TileDB-sumstat/blob/main/CONTRIBUTING.md) on how to contribute to this repo.

#### Additional link to help you get to know TileDB more

- [TileDB](https://tiledb.com/)

---
### Installation

To install this program run the following:

#### Create a virtual environment (recommended)
conda create --name tdbsumstat --file conda-linux-64.lock

conda activate tdbsumstat

#### Install the package
make install

#### Verify installation
tdbsumstat --help

#### Run with singularity

Alternatively you can use the singularity container that already pack tdbsumstat. you can find the container here: ghcr.io/bruno-ariano/tdbsumstat:1.0


#### Run with nextflow

If you are running the pipeline with nextflow change the path where the conda env is stored in conf/base.config

#### To check how to use the Nextflow pipeline for extracting and ingesting data please refer here [HERE](https://github.com/HTGenomeAnalysisUnit/TileDB-sumstat/blob/nextflow_branch/docs/README.md)

---

### Code structure

The Python package lives under `tdbsumstat/` and is organised as follows:

```
tdbsumstat/
├── main.py                  # CLI entry-point (registers ingest + export commands)
├── cli/
│   ├── ingestion.py         # `tdbsumstat ingest` CLI command
│   └── export/              # `tdbsumstat export` CLI command (package)
│       ├── __init__.py      # re-exports the `export` command
│       ├── command.py       # CLI decorator + routing to handlers
│       ├── helpers.py       # shared TileDB open/metadata helper
│       ├── snp.py           # export by SNP list
│       ├── regions.py       # export by genomic regions
│       ├── locusbreaker.py  # locusbreaker wrapper
│       ├── metadata.py      # metadata export / recompute
│       └── traits.py        # bulk trait export
└── utils/
    ├── __init__.py          # acat_optimized, compute_pheno_variance, z_to_p_via_chi2
    ├── harmonize_ingest.py  # backward-compat shim → re-exports from utils/ingest/
    ├── ingest/              # ingestion pipeline (package)
    │   ├── __init__.py      # Harmonize class + HarmonizationError
    │   ├── errors.py        # HarmonizationError exception
    │   ├── schema.py        # SchemaMixin  – TileDB array creation
    │   ├── mapping.py       # MappingMixin – column-mapping CSV parsing
    │   ├── harmonize.py     # HarmonizeMixin – data normalisation
    │   ├── qc.py            # QCMixin       – optional gwaslab QC
    │   ├── writer.py        # WriterMixin   – TileDB data writer
    │   └── metadata.py      # MetadataMixin – metadata management
    ├── locusbreaker.py      # pandas-based locusbreaker (legacy)
    ├── locusbreaker_plpl.py # Polars-based locusbreaker (used by export)
    └── update_metadata.py   # standalone metadata update utility
scripts/
    ├── create_metadata.py           # one-off metadata creation helper
    ├── fix_json.py                  # one-off JSON repair helper
    └── generate_table_cell_sumstat.py  # one-off table generation helper
```

#### Ingestion pipeline in detail

The `Harmonize` class (in `tdbsumstat/utils/ingest/`) is composed from focused mixin classes:

| Module | Mixin | Responsibility |
|--------|-------|----------------|
| `schema.py` | `SchemaMixin` | Create the TileDB sparse array schema |
| `mapping.py` | `MappingMixin` | Parse the column-mapping CSV |
| `harmonize.py` | `HarmonizeMixin` | Rename columns, handle alleles, compute p-values |
| `qc.py` | `QCMixin` | Optional gwaslab-based QC checks |
| `writer.py` | `WriterMixin` | Deduplicate and append data to TileDB |
| `metadata.py` | `MetadataMixin` | Create, merge and export metadata JSON/CSV |

A typical ingestion workflow (Python API):

```python
from tdbsumstat.utils.ingest import Harmonize
import polars as pl

h = Harmonize(
    mapping_file="mapping.csv",
    uri="my_tiledb",
    type_sumstat="qtl",   # or "gwas"
    pvar_file=None,
    type_trait="quant",   # or "binary"
    mac=None,
    maf=None,
    permuted=False,
)

# One-time setup
h.create_tiledb()
h.create_mapping()

# Per-file loop
for filepath, cell, gene in file_list:
    sumstat = pl.read_csv(filepath, separator="\t", null_values="NA")
    h.harmonize(sumstat=sumstat, cell=cell, gene=gene, n=n, pheno_var=pheno_var)
    h.ingest_data(file_path=filepath)
    h.create_metadata(file_path=filepath)

# Finalise metadata
h.merge_metadata_files()
```

---

### Running tests

Tests live in the `tests/` directory and use [pytest](https://pytest.org).

```bash
# Install dev dependencies
pip install pytest pytest-cov

# Run the full test suite
python -m pytest tests/ -v

# Run only ingestion tests
python -m pytest tests/test_ingest.py -v

# Run only export tests
python -m pytest tests/test_export.py -v

# Run with coverage report
python -m pytest tests/ --cov=tdbsumstat --cov-report=term-missing
```

Test files:

| File | What it tests |
|------|---------------|
| `tests/test_utils.py` | `acat_optimized`, `compute_pheno_variance`, `z_to_p_via_chi2` |
| `tests/test_harmonize.py` | `Harmonize` class (legacy + backward-compat) |
| `tests/test_ingest.py` | Each ingest mixin + end-to-end pipeline with example data |
| `tests/test_export.py` | All export modules + CLI command routing |

