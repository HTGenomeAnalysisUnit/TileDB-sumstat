# TileDB for single cell analysis

<img width="345" alt="Screenshot 2024-11-17 at 20 47 08" src="https://github.com/user-attachments/assets/cd7ffb0e-1d58-41a9-9b86-d29420849111">

## This is a program for the ingestion and analysis of single cell QTL data using TileDB

### Description

Welcome aboard fellow developer, this is a journey we are taking for implementing a custom TileDB for single cell QTL data.
Please have a look [here](https://github.com/HTGenomeAnalysisUnit/TileDB-SC-QTL/blob/main/CONTRIBUTING.md) on how to contribute to this repo.

## Additional link to help you get to know TileDB, Dask and Polars more

- [TileDB](https://tiledb.com/)
- [Dask](https://www.dask.org/)
- [Polars](https://pola.rs/)

---

#### Note: When you add a new feature or modify an existing one please add it to the README guidelines [HERE](https://github.com/HTGenomeAnalysisUnit/TileDB-SC-QTL/blob/main/docs/README.md)

## Installation

To install this program run the following:

### Create a virtual environment (recommended)
conda create --name scqtl --file conda-{linux, osx-arm}-64.lock

conda activate scqtl

### Install the package
make install

### Verify installation
scqtl --help
