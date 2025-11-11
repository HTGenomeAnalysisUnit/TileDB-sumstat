# TileDB for summary statistics analysis

## Nextflow usage
This pipeline consists of 2 main subcommands called ingestion and export. The ingestion is used in import new data into an existing or a new TileDB while the export can be used to query the data and export it in a txt file. This program makes also optionally use of Nextlfow for speed-up the computation.

# Flanders : Finemapping coLocalization AND plEiotRopy Solver

Flanders is a modular pipeline and toolkit for scalable **fine-mapping** and **colocalization** of genetic association signals across large scale datasets and multiple traits.
Implemented using **Nextflow** and mostly **R**, it separates computationally intensive fine-mapping from downstream colocalization to optimize reusability and performance.

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
### Example: Fine-Mapping + Colocalization
```bash
nextflow run main.nf -r 1.0    -profile [docker|singularity|conda]    --summarystats_input /path/to/input_table.tsv    --run_colocalization true    --finemap_id my_finemap_run    --coloc_id my_coloc_run    -w ./work    -resume
```

### Example: Run Only Colocalization (with existing `.h5ad`)
```bash
nextflow run Biostatistics-Unit-HT/Flanders -r 1.0    -profile [docker|singularity|conda]    --coloc_h5ad_input /path/to/finemapping_output.h5ad    --run_colocalization true    --coloc_id my_coloc_run    -w ./work    -resume
```

### Quick run with example dataset
```bash
nextflow run Biostatistics-Unit-HT/Flanders -r 1.0 -profile test,singularity -w ./work
```
</br>

## 🧠 Pipeline overview
Flanders separates the fine-mapping and colocalization process into two distinct steps:
</br>
</br>
### Step 1: Fine-mapping

#### Required Inputs
 Input | Description |
|-------|-------------|
| [GWAS summary statistics](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/Fine‐mapping-inputs#gwas-summary-statistics) | `.tsv`/`.csv` (optionally gzipped) |
| LD reference panel | PLINK-format reference panel (`.bed/.bim/.fam`) — preferably from the same sample population used in the GWAS |
| [Metadata and GWAS-specific parameters table](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/Fine‐mapping-inputs#metadata-and-gwas-specific-parameters-table) | `.tsv` file listing GWAS summary statistics paths and trait-specific parameters |
</br>

#### Steps

1. **Munging of GWAS summary statistics**
  - Format harmonization and imputation of missing information (e.g. missing allele frequency calculated from the LD reference panel)
  - Optional liftover to GRCh38
  - Optional restriction of analysis to enlisted chromosomes
    </br>⚠️ In addition to autosomes, chromosomes **X** and **Y** are also accepted.
  - Alphabetical ordering of alleles, ensuring the first one in alphabetical order is the effect allele (effect sizes and allele frequencies are flipped/inverted where needed)
  - Conversion of SNP IDs to Flanders internal coding of `"chr"CHR:POS:EA:NEA` ***where EA is the first allele in alphabetical order***
    </br>⚠️ This differs from common REF/ALT conventions and allows for robust variants matching between multiple GWAS summary statistics and LD reference panel.
</br>

#### 2. Identification of significantly associated genomic regions
  - Identifies genomic regions containing significant associated SNPs by employing `Locusbreaker`, an in-house developed algorithm which defines each association peak based on the   distance between the end of a peak and the start of the next one.

<details>
  
`Locusbreaker` first selects all SNPs below a given a p-value threshold (suggested value 1x10<sup>-6</sup>, customizable at the column `p_thresh2` of the [Metadata and GWAS-specific parameters table](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/Fine‐mapping-inputs#metadata-and-gwas-specific-parameters-table)), identifying groups of SNPs positionally close to each other.
</br>If two consecutive SNPs are closer to each other than a set distance threshold (suggested value 250kb, customizable at the column `hole` of the [Metadata and GWAS-specific parameters table](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/Fine‐mapping-inputs#metadata-and-gwas-specific-parameters-table)), they are grouped into the same locus, while if they are further apart than the distance threshold, they are used to define the boundaries between peaks.
</br>Loci with at least a significant SNPs (suggested value 5x10<sup>-8</sup>, customizable at the column `p_thresh1` of the [Metadata and GWAS-specific parameters table](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/Fine‐mapping-inputs#metadata-and-gwas-specific-parameters-table)) are retained and their boundaries are enlarged by 100kb to fully capture the shape of the association peak.
</details>
  
</br>


#### 3. Fine-Mapping with SuSiE-RSS
  - For each genomic region, finemapping is performed using [SuSiE-RSS](https://stephenslab.github.io/susieR/reference/susie_rss.html) and LD calculated from input PLINK files
</br>⚠️ Whenever possible, in sample LD is strongly recommended (especially for molecular omic phenotypes where the explained variance can be very large).
</br>⚠️ Be aware that only SNPs in common between the GWAS summary statistics and the LD reference panel are taken into account for fine-mapping, while all other SNPs are discarded (loci for which no SNP overlap is found between the GWAS summary statistics and the LD reference panel are reported in [NOT_FINEMAPPED_no_variants_from_locus_in_LD_ref.tsv](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/Fine%E2%80%90mapping-exceptions#not_finemapped_no_variants_from_locus_in_ld_reftsv)).
</br>⚠️ Be aware that **loci fully or partially overlapping the HLA region (GRCh38: chr6:28,510,120-33,480,577) are excluded from fine-mapping**. The HLA region is characterized by extremely high variant density, long-range linkage disequilibrium and complex haplotype patterns, which can bias statistical fine-mapping methods and reduce confidence in inferred causal variants.
</br>

#### 4. Saving fine-mapping results to AnnData object
  - Log approximate Bayes factors (lABFs) and metadata for the 99% credible sets are stored in an [AnnData object](https://anndata.dynverse.org/index.html) (.h5ad).
</br>
</br>

### Step 2: Colocalization analysis

#### Inputs
|Input | File description |
|------|------------------|
| [Fine-mapping AnnData](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/csAnnData-specifications)    |	An `.h5ad` file containing lABFs and metadata of credible sets (output from the fine-mapping step)  |
</br>

#### Steps

#### 1. Generation of colocalization guide table
- Lists all pairs of credible sets that share at least one SNP  (it is not possible for credible sets to colocalize without sharing at least a SNP).
  </br>⚠️ If no credible sets share at least one SNP, no colocalization is performed and an empty guide table is produced.

#### 2. Colocalization with iCOLOC
  - Performs pair-wise colocalization for pair of credible sets listed in the guide table by employing `iCOLOC`, a framework extending traditional [colocalization analysis using Bayes Factors](https://chr1swallace.github.io/coloc/reference/coloc.abf.html) by imputing lABFs of SNPs outside of credible sets to the minimum lABF value in the locus.

<details>
  
iCOLOC approach allows to:
  1. Significantly reducing storage requirements by saving in the AnnData object only exact lABF values of credible sets SNPs
  2. Enhancing colocalization accuracy compared to tradional coloc by reducing false positives due to two causal SNPs being in strong LD.
</details>

</br>
</br>


## 📁 Output
| Output Type                                     | Description                                                                                                      |
| ----------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| `gwas_and_loci_tables/*_dataset_aligned.tsv.gz` | Harmonized (and optionally lifted) GWAS summary statistics                                                       |
| `gwas_and_loci_tables/*_loci.tsv`               | Boundaries of identified association regions and GWAS summary statistics for the sentinel SNP                    |
| `finemapping_exceptions/`                       | [Multiple tables reporting information about loci that were not fine-mapped with the standard procedure or at all](https://github.com/Biostatistics-Unit-HT/Flanders/wiki/Fine%E2%80%90mapping-exceptions) |
| `finemapping/*_susie_finemap.rds` _(optional)_  | Individual RDS files for each fine-mapped locus                                                                  |
| `anndata/*.h5ad`                                | AnnData object with lABFs, CS metadata and SNP annotations resulting from fine-mapping                           |
| `coloc/coloc_guide_table.csv`                   | Colocalization analysis guide table, listing all colocalization tests performed                                  |
| `coloc/*_colocalization.table.*.tsv`            | Colocalization analysis results (all, filterd by PPH4 threshold and filtered by PPH3 threshold)                  |
</br>
</br>


## 👩‍🔬 Credits
Developed by the Biostatistics and Genome Analysis Units at [Human Technopole](https://humantechnopole.it/en/)<br>
-  [Arianna Landini](mailto:arianna.landini@fht.org)<br>
-  [Sodbo Sharapov](mailto:sodbo.sharapov@fht.org)<br>
-  [Edoardo Giacopuzzi](mailto:edoardo.giacopuzzi@fht.org)<br>
-  [Bruno Ariano](mailto:bruno.ariano@fht.org)<br>
-  [Nicola Pirastu](mailto:nicola.pirastu@fht.org)<br>

## TileDB advanced usage

```
python main.py

Usage: main.py [OPTIONS] COMMAND [ARGS]...

  Single cell analysis with TileDB

Dask cluster options:
  --workers INTEGER  Number of workers to use in the Dask cluster
  --cores_w INTEGER  Number of core to use per single worker
  --memory_w TEXT    Number of GB to give as memory for each worker
  --local_cluster    Use this option to parallelize within a node (Use for Sanger cluster)

Other options:
  --help             Show this message and exit.

Commands:
  ingestion  Ingest single cell QTL with TileDB
  export     Query the TileDb and export the data in csv format
```
</br>

test
## Ingestion
To ingest new data into a TileDB use the ```ingestion``` option of the program. To check all the possible option that can be given below.

``` python main.py ingestion --help```

<details>
  <summary>Help message</summary>
  
  ```
  Ingest single cell QTL with TileDB

  Essential parameters:
  --uri TEXT         Where to store the TileDB
  --list_files TEXT  List of the files to ingest. Each files needs to be in tsv format (gz or not)
  --cell_type TEXT   The celltype to ingest

  Optional parameters:
  --batch_size INT  The number of rows form a file to ingest at once(Don't touch this parameter unless you know how ot will affect both the memory of the process and the TileDB itself).
  --chunk_size INT  The number of rows to ingest at once
  Other options:
  --help             Show this message and exit.
  ```
</details>
</br>

### Query and Export
To query and export the TileDB you need to use the command exportlike shown below. Also below are the options you can see the type of query you can use and how to export the data in a csv file

``` python main.py export --help```

<details>
  <summary>Help message</summary>

```
Usage: main.py export [OPTIONS]

  Qeury TileDB database and export data.

Options for querying the TileDB:
  --uri   TEXT            Where to data to be created or queried is stored
  --schema                Print the schema of a tiledb
  --chrom TEXT            Cell to interrogate
  --cell  TEXT            Cell to interrogate
  --gene  TEXT            Gene to interrogate
  --snp   TEXT            List of SNPs to interrogate taken from a txt file.
                          Please check the example file in test_data/dummy_SNPs.txt for details on the format of this file
  

Options for general filters into TileDB:
  
  --maf   FLOAT   

Options for Locusbreaker:
  --locusbreaker        Option to run locusbreaker
  --pvalue-sig    FLOAT      P-value threshold to use for filtering the data (default: 5e-8)
  --pvalue-limit  FLOAT      P-value threshold for loci borders (default: 5e-6)
  --hole-size     INTEGER    Minimum pair-base distance between SNPs in different loci (default: 250000)
  --maf           FLOAT      The MAF to filter the TILEDB before locusbreaker (default: 0.01)
  --table         TEXT       Path of the table containing the traits to interrogate

Options for output:
  --out           TEXT       Output path with file name where results will be stored

                        
Other options:
  --help              Show this message and exit.
```
</details>
</br>

## Dry example run 

To give an idea on how this program works we created a series of script and dummy data for you to test the program. Please keep in mind that for running these example you need at least 15 GB of memory on your computer and you need also the conda environment activated

</br>

### Step 1 ingest your data in a TileDB

To test how to ingest data you can use the below command which use all the options described above:

```
python main.py ingestion --uri dummy_tiledb_CD14 --list_files test_data/list_files.txt
```
</br>

### Step 2 investigate the schema of the TileDB

To get details on how the TileDB you created is made run the following command:

```
python main.py export --uri dummy_tiledb_CD14 --schema
```
<details>
  <summary>Schema of the file</summary>

```
ArraySchema(
  domain=Domain(*[
    Dim(name='cell_type', domain=('', ''), tile=None, dtype='|S0', var=True, filters=FilterList([LZ4Filter(level=5), ])),
    Dim(name='gene', domain=('', ''), tile=None, dtype='|S0', var=True, filters=FilterList([LZ4Filter(level=5), ])),
    Dim(name='position', domain=(1, 3000000000), tile=3000000000, dtype='uint32', filters=FilterList([LZ4Filter(level=5), ])),
  ]),
  attrs=[
    Attr(name='SNP', dtype='ascii', var=True, nullable=False, enum_label=None, filters=FilterList([LZ4Filter(level=5), ])),
    Attr(name='af', dtype='float32', var=False, nullable=False, enum_label=None, filters=FilterList([LZ4Filter(level=5), ])),
    Attr(name='beta', dtype='float32', var=False, nullable=False, enum_label=None, filters=FilterList([LZ4Filter(level=5), ])),
    Attr(name='se', dtype='float32', var=False, nullable=False, enum_label=None, filters=FilterList([LZ4Filter(level=5), ])),
    Attr(name='allele0', dtype='ascii', var=True, nullable=False, enum_label=None, filters=FilterList([LZ4Filter(level=5), ])),
    Attr(name='allele1', dtype='ascii', var=True, nullable=False, enum_label=None, filters=FilterList([LZ4Filter(level=5), ])),
    Attr(name='p-value', dtype='float64', var=False, nullable=False, enum_label=None, filters=FilterList([LZ4Filter(level=5), ])),
  ],
  cell_order='row-major',
  tile_order='row-major',
  capacity=10000,
  sparse=True,
  allows_duplicates=True,
)
```
</details>
This will print the schema of the TileDB created.
</br>

### Step 3 query and export data by a list of SNPs

Here we are showing how to query the TileDB created above for a series of SNPs within cell types

```
python main.py export --uri dummy_tiledb_CD14 --snp test_data/dummy_SNPs.txt --cell_types test_data/dummy_cell_list.txt --output_path test_out_snp_CD14
```
</br>

### Step 4 query and export data by a list of genes and cell types for all the positions

```
python main.py export --uri dummy_tiledb_CD14 --genes test_data/dummy_gene_list.txt --cell_types test_data/dummy_cell_list.txt --output_path test_out_genes_CD14
```
</br>

### Step 5 run locus-breaker

```
python main.py export --uri dummy_tiledb_CD14 --genes test_data/dummy_gene_list.txt --cell_types test_data/dummy_cell_list.txt --output_path test_out_genes_CD14
```
This will create 2 files in txt format. One containing the the limit of each clumped region and another file with all the positions included in each independent region defined by significant SNPs.

### Step 6 run locus-breaker with Dask

To run both the ingestion and locusbreaker on parallel workers using Dask you can use the --workers parameters to define how many you want to use like described in the command below.
```
python main.py --workers 6 --memory_w 15 export --uri dummy_tiledb_CD14 --genes test_data/dummy_gene_list.txt --cell_types test_data/dummy_cell_list.txt --output_path test_out_genes_CD14
```