# TileDB-sumstat

TileDB-sumstat is a Nextflow + Python toolkit for scalable ingestion and export of genetic association summary statistics (GWAS and single‑cell QTL). It uses TileDB as the underlying storage engine and provides pipelines and utilities to import, query and export summary statistics at scale.

## Table of Contents

- [Requirements](#requirements)
- [Pipeline Overview](#pipeline-overview)
- [Usage with Nextflow](#usage-with-nextflow)
  - [Ingestion](#ingestion)
  - [Export](#export)
    - [SNP-based Export](#snp-based-export)
    - [Region-based Export](#region-based-export)
    - [Traits Export](#traits-export)
    - [Locusbreaker](#locusbreaker)
    - [Metadata Export](#metadata-export)
    - [Recompute Metadata](#recompute-metadata)
- [Nextflow Modules](#nextflow-modules)
- [Testing the Nextflow Pipeline](#testing-the-nextflow-pipeline)
  - [Stub Tests (CI / no data required)](#stub-tests-ci--no-data-required)
  - [Integration Test](#integration-test)
- [Test Data](#test-data)
- [Support and Contribution](#support-and-contribution)


---

## Requirements

Before running the pipeline, make sure you have:

- **Nextflow** (recommended v24.04+)
- **Python** (recommended 3.8+)
- **Conda** (the pipeline includes conda profiles for environment management)
- **Singularity** (the pipeline includes singularity profiles for environment management)



See the repository main README for detailed environment setup and instructions to create the conda environments or pull the singularity image used by the pipeline.

---

## Pipeline Overview

TileDB-sumstat implements two main workflows executed in separate steps:

1. **Ingestion** — Import summary statistics files into a TileDB array
2. **Export** — Query the TileDB array and export results by:
   - SNP
   - Genomic region
   - Entire summary statistics (trait-based)
   - Clumping using the "Locusbreaker" algorithm

The program can be used with Nextflow or as standalone Python utilities.

---

## Usage with Nextflow

Examples below use main.nf with placeholder values — replace them with your actual paths and parameters.

### Ingestion

Import summary statistics files into a TileDB array.

#### Required Files

- **Mapping File** (.csv): Maps input GWAS columns to standard TileDB-sumstat columns. See `example_data/mapping_file_test.csv` for format. First column: original GWAS names, second column: converted names.
- **Data Table** (`example_data/example_data_table.csv`): Lists files to ingest. Must contain **absolute file paths** in the `FILE` column. Optional columns:
  - `N`: Sample size
  - `N_CASES` / `N_CONTROLS`: Case/control sizes for binary GWAS
  - `CELL`: Cell type (single-cell QTL)
  - `GENE`: Gene name (single-cell QTL)
  - `PHENO_VAR`: Phenotypic variance (single-cell QTL)
  - `TRAIT`: Trait name (GWAS)

> **Note on file paths:** The `FILE` column in the data table must contain absolute paths that are accessible from the compute nodes running the pipeline. The file `example_data/example_data_table_test.csv` shows the column format; for production runs, replace the paths with absolute paths on your cluster.

#### Parameters

**Required:**
- `--ingestion` (flag to enable ingestion)
- `--file_path_ingestion` (path to data table file)
- `--mapping_file` (path to column mapping file)
- `--type_sumstat` (either `gwas` or `qtl`)
- `--tiledb_name` (name for the TileDB array)

**Optional:**
- `--qc` (enable QC processing via gwaslab)
- `--ingestion_chunk_files` (number of files to ingest per batch, default: 4)
- `--maf` (minimum allele frequency filter, default: 0)
- `--mac` (minimum allele count filter, default: 0)
- `--permuted` (compute SE from permuted p-value)
- `--pvar_file` (pvar file to align alleles)
- `--outdir` (output directory, default: `./results`)

#### Example
```bash
nextflow run main.nf \
  -profile singularity \
  --ingestion \
  --file_path_ingestion example_data/example_data_table.csv \
  --mapping_file example_data/mapping_file_test.csv \
  --type_sumstat qtl \
  --tiledb_name my_tiledb \
  --ingestion_chunk_files 4
```

Using a pre-defined test profile:
```bash
nextflow run main.nf -profile test_ingest,singularity
```


### Export

Query and export data from the TileDB array.

#### Common Parameters

**Required:**
- `--export` (flag to enable export)
- `--uri_path` (path to TileDB array)
- `--out` (output file prefix)
- `--type_sumstat` (type of summary statistics: `gwas` or `qtl`)

**Optional:**
- `--attrs` (attributes to export, e.g., `"BETA,SE,EAF"`, default: `"P,SNPID,EAF,BETA,SE"`)

#### SNP-based Export

Extract specific SNP positions.

**Required:**
- `--snp` (path to SNP list CSV, columns: `CHR`, `POS`, `TRAIT`)

Example files: `example_data/snp_list_sc.csv`

```bash
nextflow run main.nf \
  -profile singularity \
  --export \
  --snp example_data/snp_list_sc.csv \
  --tiledb_path /path/to/tiledb \
  --uri_path /path/to/tiledb \
  --attrs "BETA,SE,P" \
  --out results/snp_output \
  --type_sumstat qtl
```

Quick test (stub – no data needed):
```bash
nextflow run main.nf -profile test_export_snp -stub
```

#### Region-based Export

Extract genomic intervals.

**Required:**
- `--regions` (path to regions CSV, columns: `CHR`, `START`, `END`, `TRAIT`)

Example file: `example_data/region_list_sc.csv`

```bash
nextflow run main.nf \
  -profile singularity \
  --export \
  --regions example_data/region_list_sc.csv \
  --uri_path /path/to/tiledb \
  --attrs "BETA,SE,P" \
  --out results/regions_output \
  --type_sumstat qtl
```

Quick test:
```bash
nextflow run main.nf -profile test_export_regions -stub
```

#### Traits Export

Export complete summary statistics for a list of traits or cell-type/gene combinations.

**Required:**
- `--export_traits` (flag)
- `--list_traits` (path to CSV with `TRAIT` column; for QTL: `CELL:GENE` format)

Example file: `example_data/trait_list_test.csv`

```bash
nextflow run main.nf \
  -profile singularity \
  --export \
  --export_traits \
  --list_traits example_data/trait_list_test.csv \
  --uri_path /path/to/tiledb \
  --out results/traits_output \
  --type_sumstat qtl
```

Quick test:
```bash
nextflow run main.nf -profile test_export_traits -stub
```

#### Locusbreaker

Identify genomic loci with significant associations.

**Required:**
- `--locusbreaker` (flag)
- `--table_lb` (path to traits table, columns: `CHR`, `TRAIT`, `SIG`, `LIM`)

**Optional:**
- `--maf_lb` (MAF filter, default: 0.001)
- `--locus_max_size_lb` (maximum locus size in bp, default: 3 Mb)
- `--hole_lb` (maximum gap within loci in bp, default: 250 kb)
- `--cis_trans_lb` (for QTLs: `cis` or `trans`, default: `cis`)

Example:
```bash
nextflow run main.nf \
  -profile singularity \
  --locusbreaker \
  --table_lb example_data/locusbreaker_test_table_sc.csv \
  --uri_path /path/to/tiledb \
  --out results/lb_output \
  --type_sumstat qtl
```

Quick test:
```bash
nextflow run main.nf -profile test_export_lb -stub
```

#### Metadata Export

Export the merged metadata stored in the TileDB array to CSV.

```bash
tdbsumstat export --export-meta --uri-path /path/to/tiledb --type-sumstat qtl --out metadata_out
```

#### Recompute Metadata

Recompute per-trait/cell metadata statistics after applying a MAC filter without modifying the TileDB data.

**Required:**
- `--recompute_meta` (flag)
- `--list_traits` (CSV with `CELL` and `CHR` columns for QTL, or `TRAIT` and `CHR` for GWAS)
- `--uri_path`
- `--mac` (minimum allele count threshold)

Example file: `example_data/recompute_meta_test_table.csv`

```bash
nextflow run main.nf \
  -profile singularity \
  --recompute_meta \
  --list_traits example_data/recompute_meta_test_table.csv \
  --uri_path /path/to/tiledb \
  --mac 10 \
  --out results/meta_out \
  --type_sumstat qtl
```

Quick test:
```bash
nextflow run main.nf -profile test_recompute_meta -stub
```

---

## Nextflow Modules

The pipeline is composed of the following Nextflow process modules under `modules/`:

| Module | Process | Description |
|--------|---------|-------------|
| `create_tiledb/` | `CREATE_TILEDB` | Creates the TileDB sparse array schema |
| `ingestion/` | `INGEST_DATA` | Harmonises and ingests one file-list chunk into TileDB |
| `merge_metadata/` | `MERGE_METADATA` | Collects per-chunk metadata JSON files and stores merged metadata in TileDB |
| `snp/` | `EXPORT_SNP` | Exports data for a list of SNP positions |
| `regions/` | `EXPORT_REGIONS` | Exports data for a list of genomic intervals |
| `traits/` | `TRAITS` | Exports complete summary statistics for a list of traits |
| `locusbreaker/` | `EXPORT_LOCUSBREAKER` | Runs the Locusbreaker algorithm and exports locus/segment tables |
| `recompute_meta/` | `RECOMPUTE_META` | Recomputes metadata with a MAC filter |

Each module has a `stub` block so the workflow DAG can be validated without executing the actual commands (see [Testing](#testing-the-nextflow-pipeline)).

---

## Testing the Nextflow Pipeline

### Stub Tests (CI / no data required)

Stub tests validate the workflow DAG structure (channels, processes, outputs) without running the actual `tdbsumstat` commands. They are the recommended way to test the pipeline in CI or when data is not available.

All test profiles are defined in `conf/test.config` and registered in `nextflow.config`.

```bash
# Validate ingestion workflow
nextflow run main.nf -profile test_ingest -stub

# Validate export workflows
nextflow run main.nf -profile test_export_snp -stub
nextflow run main.nf -profile test_export_regions -stub
nextflow run main.nf -profile test_export_traits -stub
nextflow run main.nf -profile test_export_lb -stub
nextflow run main.nf -profile test_recompute_meta -stub
```

These stubs are also run automatically on every push and pull request via the [CI workflow](.github/workflows/ci.yml).

### Integration Test

To run a full integration test that actually executes the ingestion and export commands, you need:
1. `tdbsumstat` installed (`pip install -e .`)
2. Nextflow installed

```bash
# Build a CI-friendly data table with absolute paths
{
  echo "FILE,CELL,GENE,PHENO_VAR,N"
  echo "$(pwd)/example_data/dummy_out_ENSG0000010000.tsv.gz,Tgd,ENSG0000010000,1.5,4000"
  echo "$(pwd)/example_data/dummy_out_ENSG0000010001.tsv.gz,Tgd,ENSG0000010001,1.5,4000"
} > /tmp/example_data_table_ci.csv

# Run ingestion (no container – uses local tdbsumstat)
nextflow run main.nf \
  --ingestion true \
  --file_path_ingestion /tmp/example_data_table_ci.csv \
  --mapping_file "$(pwd)/example_data/mapping_file_test.csv" \
  --type_sumstat qtl \
  --tiledb_name test_ci \
  --ingestion_chunk_files 2 \
  --maf 0 \
  --mac 0 \
  --outdir ./results_ci \
  -process.container null \
  -ansi-log false
```

---

## Test Data

| File | Description |
|------|-------------|
| `example_data/dummy_out_ENSG0000010000.tsv.gz` | Example QTL summary statistics (gene ENSG0000010000) |
| `example_data/dummy_out_ENSG0000010001.tsv.gz` | Example QTL summary statistics (gene ENSG0000010001) |
| `example_data/example_data_table.csv` | Ingestion data table (absolute paths, for cluster use) |
| `example_data/example_data_table_test.csv` | Ingestion data table (filenames only, documents column format) |
| `example_data/mapping_file_test.csv` | Column mapping file for example data |
| `example_data/snp_list_sc.csv` | SNP list for SNP-based export tests |
| `example_data/region_list_sc.csv` | Region list for region-based export tests |
| `example_data/locusbreaker_test_table_sc.csv` | Traits table for Locusbreaker tests |
| `example_data/trait_list_test.csv` | Trait list (`TRAIT` column) for traits export tests |
| `example_data/recompute_meta_test_table.csv` | Cell/CHR table for recompute metadata tests |

---

## Support and Contribution

If you encounter issues, please open a GitHub issue with a reproducible example and error messages.

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request against the main branch
4. Follow the repository's contribution guidelines
