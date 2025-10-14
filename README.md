# TileDB for summary statistics analysis

<img width="345" alt="Screenshot 2024-11-17 at 20 47 08" src="https://github.com/user-attachments/assets/cd7ffb0e-1d58-41a9-9b86-d29420849111">

## This is a program for the ingestion and analysis of single cell QTL data using TileDB

### Description

Welcome to TileDB-sumstat. This is a custom TileDB for singel cell QTL and GWAS data.
Please have a look [here](https://github.com/HTGenomeAnalysisUnit/TileDB-sumstat/blob/main/CONTRIBUTING.md) on how to contribute to this repo.

## Additional link to help you get to know TileDB, Dask and Polars more

- [TileDB](https://tiledb.com/)
- [Dask](https://www.dask.org/)
- [Polars](https://pola.rs/)

---

#### Note: When you add a new feature or modify an existing one please add it to the README guidelines [HERE](https://github.com/HTGenomeAnalysisUnit/TileDB-sumstat/blob/main/docs/README.md)

## Installation

To install this program run the following:

### Create a virtual environment (recommended)
conda create --name scqtl --file conda-linux-64.lock

conda activate scqtl

### Install the package
make install

### Verify installation
scqtl --help
