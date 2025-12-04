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
    - [Locusbreaker](#locusbreaker)
- [Test Data](#test-data)
- [Support and Contribution](#support-and-contribution)


---

## Requirements

Before running the pipeline, make sure you have:

- **Nextflow** (recommended v24.04+)
- **Python** (recommended 3.8+)
- **Conda** (the pipeline includes conda profiles for environment management)


See the repository main README for detailed environment setup and instructions to create the conda environments used by the pipeline.

---

## Pipeline Overview

TileDB-sumstat implements two main workflows executed in separate steps:

1. **Ingestion** — Import summary statistics files into a TileDB array
2. **Export** — Query the TileDB array and export results by:
   - SNP
   - Genomic region
   - Entire summary statistics
   - Clumping using the "Locusbreaker" algorithm

The program can be used with Nextflow or as standalone Python utilities.

---

## Usage with Nextflow

Examples below use main.nf with placeholder values — replace them with your actual paths and parameters.

### Ingestion

Import summary statistics files into a TileDB array.

#### Required Files

- **Mapping File** (.csv): Maps input GWAS columns to standard TileDB-sumstat columns. See `example_data/mapping_file_test` for format. First column: original GWAS names, second column: converted names.
- **Data Table**: Lists files to ingest. See `example_data/example_data_table.csv`. Must include:
  - Path to GWAS files
  - Optional: Additional columns can be added to this file if they are not present in the summary statistics:
    - `N`: Sample size for specific summary statistics
    - `N_CASES`: Sample size for cases specific to a summary statistics
    - `N_CONTROLS`: Sample size for controls specific to a summary statistics
    - `CELL`: Cell name for single QTL studies
    - `GENE`: GENE name for single QTL studies
    - `PHENO_VAR`: Phenotypic variance for the specific trait (only for QTLs single cell)
    - `TRAIT`: Trait name for GWAS studies

#### Parameters

**Required:**
- `--ingest` (flag to enable ingestion)
- `--file_path_ingestion` (path to data table file)
- `--type_sumstat` (either `gwas` or `qtl`)

**Optional:**
- `--qc` (enable QC processing)
- `--ingestion_chunk_files` (number of files to ingest simultaneously, default: 4)

#### Example
```bash
nextflow run main.nf -profile conda --ingest --file_path_ingestion example_data/example_data_table.csv  --mapping_file example_data/mapping_file_test.csv --type_sumstat qtl --ingestion_chunk_files 4
```

Using Nextflow profile you can also use

```bash
nextflow run main.nf -profile test_ingest,conda
```


### Export

Query and export data from the TileDB array.

#### Common Parameters

**Required:**:
- --export (flag to enable export)
- --tiledb_path (path to TileDB array)
- --out (output file prefix)
- --type_sumstat (type of summary statistics, either gwas or qtl for single cell)

**Optional**:
- --attrs (attributes to export, e.g., "BETA,SE,PVAL,EAF,A1,A2")

#### SNP-based Export

Extract specific SNP positions.

**Required**
- --snp (path to SNP list file)

Example Files:
- GWAS: example_data/snp_list_gwas.csv
- Single-cell: example_data/snp_list_sc.csv

#### Example:
```bash
nextflow run main.nf -profile conda --export --tiledb_path /path/to/tiledb --snp /path/to/snp_list.csv --attrs "BETA,SE,PVAL,EAF,A1,A2" --out /path/to/output_prefix --type_sumstat gwas
```
#### Region-based Export

Extract genomic regions using BED format.

**Required:**
- --table-regions (path to regions file)

Example:
```bash
nextflow run main.nf -profile conda --export --tiledb_path /path/to/tiledb --table-regions /path/to/regions_table.csv --attrs "BETA,SE,PVAL" --out /path/to/output_prefix --type_sumstat gwas
```

Quick Test:
```bash
nextflow run main.nf -profile test_export_lb,conda
```

#### Locusbreaker

Identify genomic loci with significant associations and export locus-centric results.

**Required:**
- --table-lb (path to traits table, see example_data/locusbreaker_test_table.csv)

**Optional**
- --maf-lb (minor allele frequency filter)
- --locus-max-size (maximum locus size in base pairs)
- --hole-lb (maximum gap size within loci in base pairs)
- --cis-trans (for QTLs: filter by cis or trans)

Locusbreaker Algorithm:
1. Select SNPs below p-value threshold (suggested: 1e-6)
2. Group consecutive SNPs within distance threshold (suggested: 250 kb)
3. Retain groups containing at least one genome-wide significant SNP (suggested: 5e-8)
4. Expand locus boundaries by margin (e.g., +100 kb)
5. Apply additional filters (MAF, locus size, cis/trans)

---

#### Metadata extraction

Metadata are structured as json in tiledb. To extract them you can use the following command:

```bash
tdbsumstat export metadata 
```

tiledb_meta
{'traits': [], 'CELL': ['Tgd'], 'Tgd': {'20': {'ENSG0000010000': {'ACAT': 0.0, 'N': 4000.0, 'PHENO_VAR': 1.0}, 'ENSG0000010001': {'ACAT': 0.0, 'N': 4000.0, 'PHENO_VAR': 1.0}}}}
## Test Data

- example_data/snp_list.csv - Example SNP list for SNP-based export
- example_data/example_data_table.csv - Example table for regions and ingestion
- example_data/mapping_file_test - Example mapping file for ingestion
- example_data/locusbreaker_test_table.csv - Example table for Locusbreaker

---

## Support and Contribution

If you encounter issues, please open a GitHub issue with a reproducible example and error messages.

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request against the main branch
4. Follow the repository's contribution guidelines

---