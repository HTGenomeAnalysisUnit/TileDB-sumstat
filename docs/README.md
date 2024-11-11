# TileDB for single cell analysis

## Overall usage
This program works using. python on the Sanger cluster using the tiledb conda environment. This program consists of 2 main subcommands called ingestion and export. The ingestion is used in import new data into an axisting or a new TileDB while the export can be used to query the data and export in different format. This program makes also optionally use of Dask for accelerate the computation. To get a menu of the main parameters and commands available run a script like below.

```
$user12 python main.py

Usage: main.py [OPTIONS] COMMAND [ARGS]...

  Single cell analysis with TileDB

Dask cluster options:
  --workers INTEGER  Number of workers to use in the Dask cluster
  --cores_w INTEGER  Number of core to use per single worker
  --memory_w TEXT    Number of GB to give as memory for each worker

Other options:
  --help             Show this message and exit.

Commands:
  ingestion  Ingest single cell QTL with TileDB
  export  Query the TileDb and export the data in csv format
```
</br>


## Ingestion
To ingest new data into a TileDB use the ingestion option of the program. To check all the possible option run

``` python main.py ingestion --help```

<details>
  <summary>Help message</summary>
  
  ```
  Ingest single cell QTL with TileDB

  Essential parameters:
  --uri TEXT         Where to store the TileDB
  --input TEXT       A tsv file in gzip format to ingest
  --cell_type TEXT   The celltype to ingest

  Optional parameters:
  --chunk_size TEXT  The number of rows to ingest at once

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
  --uri TEXT          Where to data to be created or queried is stored
  --schema            Print the schema of a tiledb
  --cell_types TEXT   List of cells to interrogate taken from a txt file
  --genes TEXT        List of genes taken from a txt file
  --snp TEXT          List of SNPs to interrogate taken from a txt file. Please
                      check README for details on the format of this file
  --output_path TEXT  Output path with file name where results will be stored

Other options:
  --help              Show this message and exit.
```
</details>
</br>


## Dry example run 

To give an idea on how this program works we created a series of script and dummy data for you to test the program. Please keep in mind that for running these example you need at least 15 GB of memory on your computer

</br>

### Step 1 ingest your data in a TileDB

To test how to ingest data you can use the below command which use all the options described above:

```
python main.py ingestion --uri dummy_tiledb_CD14 --input test_data/dummy_out.tsv.gz --cell_type CD14
```
This should create a TileDB in a few minutes (depending on how fast your file system is in terms of IO).

</br>

### Step 2 investigate the schema of the TileDB

To get details on how the TileDB you created is made run the following command:

```
python main.py export --uri dummy_tiledb_CD14 --schema
```
</br>

### Step 3 query and export data by a list of SNPs

Here we are showing how to query the TileDB created above for a series of SNPs within cell types

```
python main.py export --uri dummy_tiledb_CD14 --snp test_data/dummy_SNPs.txt --cell_types test_data/dummy_cell_list.txt --output_path test_out_snp_CD14.csv
```
</br>

### Step 4 query and export data by a list of genes and cell types for all the positions

```
python main.py export --uri dummy_tiledb_CD14 --genes test_data/dummy_gene_list.txt --cell_types test_data/dummy_cell_list.txt --output_path test_out_genes_CD14.csv
```