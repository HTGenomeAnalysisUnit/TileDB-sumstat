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
    - [Trait-based Export](#trait-based-export)
    - [Locusbreaker](#locusbreaker)
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
   - Trait (the entire summary statistics of a list of traits or cell type/gene pairs)
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
nextflow run main.nf -profile singularity --ingest --file_path_ingestion example_data/example_data_table.csv  --mapping_file example_data/mapping_file_test.csv --type_sumstat qtl --ingestion_chunk_files 4
```

Using Nextflow profile you can also use

```bash
nextflow run main.nf -profile test_ingest,singularity
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
- --attrs (attributes to export, e.g., "BETA,SE,P,EAF,SNPID")

> Note: the SNP export reads the array location from `--tiledb_path`, while the trait export
> and Locusbreaker read it from `--uri_path`. The trait export and Locusbreaker also run
> without `--export`, which only gates the SNP export.

#### SNP-based Export

Extract specific SNP positions.

**Required**
- --snp (path to SNP list file)

Example Files:
- GWAS: example_data/snp_list_gwas.csv
- Single-cell: example_data/snp_list_sc.csv

#### Example:
```bash
nextflow run main.nf -profile singularity --export --tiledb_path /path/to/tiledb --snp /path/to/snp_list.csv --attrs "BETA,SE,PVAL,EAF,A1,A2" --out /path/to/output_prefix --type_sumstat gwas
```
#### Region-based Export

Extract genomic regions using BED format.

**Required:**
- --table-regions (path to regions file)

Example:
```bash
nextflow run main.nf -profile singularity --export --tiledb_path /path/to/tiledb --table-regions /path/to/regions_table.csv --attrs "BETA,SE,PVAL" --out /path/to/output_prefix --type_sumstat gwas
```

Quick Test:
```bash
nextflow run main.nf -profile test_export_lb,singularity
```

#### Trait-based Export

Export the complete summary statistics of a list of traits (GWAS) or of a list of cell
type/gene pairs (single-cell QTL). Every association of every trait in the list is exported,
so use the SNP or region export instead when only part of a signal is needed.

**Required:**
- --export_traits (flag to enable the trait export)
- --uri_path (path to the TileDB array)
- --list_traits (path to the trait list, the file name must end in `.csv`)
- --type_sumstat (either gwas or qtl for single cell)

**Optional:**
- --attrs (columns to export, default "P,SNPID,EAF,BETA,SE")
- --out (file name prefix of the exported tables, default "out")
- --tiledb_batch_size (number of traits exported per parallel job, default 4)
- --outdir (directory where results are published, default "./results")

The trait export does not need the `--export` flag, and it reads the array from `--uri_path`
rather than `--tiledb_path`. Any `--type_sumstat` value other than `gwas` is handled as a QTL
array, so `qtl` and `eqtl` behave identically.

##### Trait list

A CSV file with a `TRAIT` column; any other column is ignored. For GWAS arrays `TRAIT` is the
trait name used at ingestion, for single-cell QTL arrays it is `CELL~GENE` (the cell type and
the gene separated by a tilde):

```csv
TRAIT
Tgd~ENSG0000010000
Tgd~ENSG0000010001
```

> Careful: the tilde is specific to the trait export. The SNP export, the region export and
> Locusbreaker still expect their cell type and gene to be separated by a colon (`CELL:GENE`),
> so the same list file cannot be reused across them as is.

The list is split into chunks of `--tiledb_batch_size` rows and one job is submitted per chunk,
so a list of thousands of genes is exported in parallel.

##### Exported columns

`--attrs` accepts both dimensions and attributes of the array:

- dimensions: `CHR`, `CELL`, `GENE`, `POS` (QTL) or `CHR`, `TRAIT`, `POS` (GWAS)
- attributes: `SNPID`, `RSID`, `EAF`, `BETA`, `SE`, `P` and, for QTL arrays, `DIST`

Alleles are not stored as separate columns, they are part of `SNPID` (`CHR:POS:A1:A2`).
Asking for a name that does not exist in the array stops the job with the list of valid names.

Each job writes one CSV per batch, `<out>_<batch index>.csv`, published under
`${outdir}/gwas_and_loci_tables/`. The files have **no header line**: the columns are the
dimensions that were not listed in `--attrs`, followed by the columns of `--attrs` in the order
they were given. Note that `--out` is a file name prefix and not a directory, the output
location is controlled by `--outdir`.

#### Example:
```bash
PIPELINE="/path/to/TileDB-sumstat/main.nf"

TILEDB_PATH="/path/to/tiledb_array"
TRAIT_CSV="/path/to/trait_list.csv"

ATTRS="SNPID,CHR,POS,DIST,EAF,BETA,SE,P"
TYPE_SUMSTAT="qtl"
OUT_PREFIX="trait_export"
OUTDIR="./results"

nextflow run "$PIPELINE" \
  --export_traits \
  --uri_path "$TILEDB_PATH" \
  --list_traits "$TRAIT_CSV" \
  --attrs "$ATTRS" \
  --out "$OUT_PREFIX" \
  --outdir "$OUTDIR" \
  --type_sumstat "$TYPE_SUMSTAT" \
  --tiledb_batch_size 100 \
  -profile singularity \
  -resume
```

On an HPC cluster, add your institutional configuration with `-c /path/to/your.config` and the
matching `-profile`.

Quick Test:
```bash
nextflow run main.nf -profile test_export_traits,singularity
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