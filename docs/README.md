# TileDB for summary statistics analysis

## Nextflow usage
This pipeline consists of 2 main subcommands called ingestion and export. The ingestion is used in import new data into an existing or a new TileDB while the export can be used to query the data and export it in a txt file. This program makes also optionally use of Nextlfow for speed-up the computation.

# TileDB-sumstat

TileDB-sumstat is a pipeline and toolkit for scalable **ingestion** and **export** of genetic association signals across large scale datasets for both GWAS and QTL(single cell) data.
Implemented using **Nextflow** **Python** and **TileDB**

---

## ✅ Requirements
Before running the pipeline, ensure you have the following installed:

- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) (v24.04+)
- For environment management, one of:
  - [Docker](https://www.docker.com/)
  - [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html)
  - [Conda](https://docs.conda.io/en/latest/)
</br>


## ▶️ Running the pipeline
### Example: Run SNPs extraction
```bash
nextflow run main.nf -profile [docker|singularity|conda] --export true --tiledb_path /path/to/tiledb --snp /path/to/snp --attrs attributes_list_string --out "out" --type_sumstat gwas(or qtl) --tiledb_batch_size number_of_batches_per_job

```

### Example: Run ingestion of sumstat
```bash
nextflow run main.nf -profile [docker|singularity|conda] --ingest true --file_path_ingestion /path/to/file_ingestion    --mapping_file path_for_mapping_file  --type_sumstat gwas(or qtl) --qc false --ingestion_chunk_files 2 --tiledb_name "test"
```

### Quick run with example dataset
```bash
nextflow run main.nf -profile test_export_lb,conda 
```
</br>

## 🧠 Pipeline overview
TileDB-sumstat implement 2 methods for ingesting end querying summary statistics in 2 disctinct steps.
</br>
</br>
### Step 1: Ingestion
This part is to ingest summary in TileDB

#### Required Inputs
Parameters | Description |
  |-------|-------------|
  | [file_path_ingestion] [List of GWAS summary statistics] A `.csv` file containing the path the summary statistics to ingests
  | [mapping_file] [Mapping file] A mapping file containing the name of the column and supposed mapped ones
  | [type_sumstat] [Type of summary statistics] This can be either a gwas or a qtl(for single cell)
  | [qc] [Perform QC] This parameter perform QC on summary statistics prior ingestion
  | [ingestion_chunk_files] [Number of summary stats in parallel] This dictate how many summary statistics to ingest in parallel
  | [type_sumstaat] [gwas, qtl]

** Note**
Check the example files in the folder
</br>

### Step 2: Export
This part is to export summary in TileDB. There is in this step a series of common parameters and a series of specific one depending on which type of export you are doing (SNP/region or LocusBreaker)

#### Common parameters
Parameters | Description |
  | [export] [Activate export] A `.csv` file containing the path the summary statistics to ingests
  | [tiledb_path] [Path for TileDB] A mapping file containing the name of the column and supposed mapped ones
  | [attrs] [BETA,SE,PVAL,EAF,A1,A2] Attributes that are used in the export
  | [tiledb_batch_size] [Divide the entire set of SNPs in batches] This dictate how many summary statistics to ingest in parallel
  | [out] [Name of the output file] Depending on the batch defined a suffix is attached to this.


#### 2.1. Export SNP/Regions
- There are 2 options to export raw sumstat, either SNPs or regions. 
Parameters | Description |
  |-------|-------------|
  | [snp] [Path of SNP] This can be either a gwas or a qtl(for single cell)


#### 2.1. Export Locusbreaker
- Identifies genomic regions containing significant associated SNPs by employing `Locusbreaker`, an in-house developed algorithm which defines each association peak based on the   distance between the end of a peak and the start of the next one.

Parameters | Description |
  |-------|-------------|
  | [maf_lb] The minor allele frequency filter to apply for locus-breaker
  | [locus_max_size_lb] Max size of the LB region
  | [hole_lb] Hole that define the distribution of the peak
  | [cis_trans_lb] In case of qtl if filtering for either cis or trans (1Mb distance)
  | [table_lb] A table containing the traits to perform LB on (check example_data/locusbreaker_test_table.csv)

<details>
`Locusbreaker` first selects all SNPs below a given a p-value threshold (suggested value 1x10<sup>-6</sup>, customizable at the column `LIM` of the [table_lb] file), identifying groups of SNPs positionally close to each other.
</br>If two consecutive SNPs are closer to each other than a set distance threshold (250kb), they are grouped into the same locus, while if they are further apart than the distance threshold, they are used to define the boundaries between peaks.
</br>Loci with at least a significant SNPs (suggested value 5x10<sup>-8</sup>, customizable at the column `SIG` of the [table_lb] file). are retained and their boundaries are enlarged by 100kb to fully capture the shape of the association peak.
</details>
</br>
