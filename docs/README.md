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
```

## Ingestion
To ingest new data into a TileDB use the ingestion option of the program. To check all the possible option run
``` python main.py ingestion```
<details>
  <summary>Help message</summary>
  
  ```
  $user12 python main.py ingestion

  Ingest single cell QTL with TileDB

  Essential parameters:
  --uri TEXT         Where to store the TileDB
  --input TEXT       A tsv file to ingest
  --cell_type TEXT   The celltype to ingest

  Optional parameters:
  --chunk_size TEXT  The number of rows to ingest at once

  Other options:
  --help             Show this message and exit.
  ```
</details>


Usage example
  ```
  $user12 python main.py ingestion --uri  /Users/bruno.ariano/work/HT/sanger_work/tiledb-scqtl/test_tiledb --input test_data.txt.gz --cell_type CD14
  ```


### Query
In development.


### Dummy file creation

To create a dummy file for test like the one on Sanger you can use the script in utils/create_dummy_data.py. To run this program you can run the script like below

```
  $user12 python utils/create_dummy_data.py
  ```