# TileDB for single cell analysis

## Overall usage
This program works using. python on the Sanger cluster using the tiledb conda environment. This program consists of 2 main subcommands called ingestion and export. The ingestion is used in import new data into an existing or a new TileDB while the export can be used to query the data and export it in a txt file. This program makes also optionally use of Dask for accelerate the computation. To get a menu of the main parameters and commands available run a script like below.

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