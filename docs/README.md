# TileDB-sumstat

TileDB-sumstat is a Nextflow + Python toolkit for scalable ingestion and export of genetic association summary statistics (GWAS and single‑cell QTL). It uses TileDB as the underlying storage engine and provides pipelines and utilities to import, query and export summary statistics at scale.

---

## Requirements

Before running the pipeline, make sure you have:

- Nextflow (recommended v24.04+)
- Python (recommended 3.8+)
- TileDB (library/runtime as required by the Python TileDB package)
- Conda (the pipeline includes conda profiles for environment management)
- A working Slurm cluster (only if you use the `conda_slurm` profile) or a local environment

See the repository main README for detailed environment setup and instructions to create the conda environments used by the pipeline.

---

## Quick start

The pipeline exposes two primary workflows: ingestion and export. You can run either (or both) using Nextflow.

Examples below use `main.nf` and placeholder values — replace them with your actual paths and parameters.

Example: export SNPs
```bash
nextflow run main.nf -profile conda_slurm \
  --export true \
  --tiledb_path /path/to/tiledb \
  --snp /path/to/snp_list.csv \
  --attrs "BETA,SE,PVAL,EAF,A1,A2" \
  --out /path/to/output_prefix \
  --type_sumstat gwas \
  --tiledb_batch_size 100
```

Example: export regions
```bash
nextflow run main.nf -profile conda \
  --export true \
  --tiledb_path /path/to/tiledb \
  --table-regions /path/to/regions_table.csv \
  --attrs "BETA,SE,PVAL" \
  --out /path/to/output_prefix \
  --type_sumstat gwas \
  --tiledb_batch_size 100
```

Example: ingestion of summary statistics
```bash
nextflow run main.nf -profile conda \
  --ingest true \
  --file_path_ingestion /path/to/files_table.csv \
  --mapping_file /path/to/mapping_file.tsv \
  --type_sumstat gwas \
  --qc false \
  --ingestion_chunk_files 4
```

- A mapping file maps observed column names in your input files to the standardized column names expected by TileDB-sumstat. See `example_data/mapping_file_test` for an example.
- A table listing files to ingest is provided in `example_data/example_data_table.csv`.
- Example SNP list: `example_data/snp_list.csv`.
- Example regions table: `example_data/example_data_table.csv`.

Quick test run (uses small test datasets/configs):
```bash
nextflow run main.nf -profile test_export_lb,conda
```

---

## Pipeline overview

TileDB-sumstat implements two main workflows executed in separate steps:

1. Ingestion — import summary statistics files into a TileDB array.
2. Export — query the TileDB array and export results by SNP, region or using the "Locusbreaker" algorithm.

### Step 1 — Ingestion

Purpose: ingest summary statistics files into TileDB.

Required inputs and key parameters:

| Parameter | Description |
|---|---|
| file_path_ingestion | CSV file listing summary-statistic files to ingest (paths). |
| mapping_file | Mapping file that maps columns in each input file to the standardized column names required by TileDB-sumstat. |
| type_sumstat | Type of summary statistics: `gwas` or `qtl` (single-cell QTL). |
| qc | Whether to perform QC before ingestion (`true`/`false`). |
| ingestion_chunk_files | Number of files to ingest in parallel (controls parallelism). |

Notes:
- Check `example_data/` for sample ingestion tables and mapping files.
- Ensure the mapping file correctly maps input column names (e.g., SNP, CHR, POS, A1, A2, BETA, SE, PVAL, etc.) to expected internal names.

### Step 2 — Export

Purpose: query TileDB and export summary statistics. There are common parameters and mode-specific parameters (SNP, region, or Locusbreaker).

Common parameters:

| Parameter | Description |
|---|---|
| export | Activate export workflow (`true`/`false`). |
| tiledb_path | Path to the TileDB dataset/array. |
| attrs | Comma-separated list of attributes to export (e.g., `BETA,SE,PVAL,EAF,A1,A2`). |
| tiledb_batch_size | Batch size for processing the set of SNPs or regions (controls memory/parallelism). |
| out | Output path or prefix. Batch runs append a suffix for each batch. |

Export modes

- SNP/Regions export:
  - SNP export: provide `--snp /path/to/snp_list.csv` (list of SNP identifiers or positions).
  - Regions export: provide `--table-regions /path/to/regions_table.csv` (table with regions to query).

- Locusbreaker export:
  - Locusbreaker identifies genomic loci with significant associations and exports locus-centric results (peaks).
  - Parameters:

| Parameter | Description |
|---|---|
| maf_lb | Minor allele frequency filter applied before locus calling. |
| locus_max_size_lb | Maximum allowed size for a Locusbreaker region (e.g., in base pairs). |
| hole_lb | Maximum allowed gap (in base pairs) inside a locus before splitting (defines peak continuity). |
| cis_trans_lb | For QTLs: choose whether to filter by cis or trans (e.g., cis = within 1Mb). |
| table_lb | Table listing traits / datasets to perform locus-breaking on (see `example_data/locusbreaker_test_table.csv`). |

Locusbreaker algorithm (brief)
- Select SNPs below a given p-value threshold (suggested LIM = 1e-6). SNPs are grouped if consecutive SNPs are closer than a distance threshold (suggested 250 kb).
- Groups (putative loci) are retained if they contain at least one genome-wide significant SNP (suggested SIG = 5e-8).
- Locus boundaries can be expanded by a margin (e.g., +100 kb) to ensure full coverage of the association signal.
- Additional filters (MAF, maximum locus size, cis/trans selection for QTL) are applied during locus selection.

---

## Examples and test data

- example_data/snp_list.csv — example SNP list for SNP-based export.
- example_data/example_data_table.csv — example table used for regions and ingestion.
- example_data/mapping_file_test — example mapping file for ingestion.
- example_data/locusbreaker_test_table.csv — example table for running Locusbreaker.

---

## Support and contribution

If you find issues, please open an issue in the repository with a reproducible example and any error messages. Contributions are welcome — please open a pull request against the `main` branch and follow the repository contribution guidelines.

---

If you'd like, I can open a pull request that replaces the current docs/README.md with this cleaned and corrected version. Would you like me to create that PR? 
```