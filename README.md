# TileDB for summary statistics analysis

<img width="345" alt="Screenshot 2024-11-17 at 20 47 08" src="https://github.com/user-attachments/assets/cd7ffb0e-1d58-41a9-9b86-d29420849111">

### This is a program for the and analysis of single cell QTL and GWAS data using Nextflow and TileDB

#### Description

Welcome to TileDB-sumstat. This is a custom TileDB for singel cell QTL and GWAS data.
Please have a look [here](https://github.com/HTGenomeAnalysisUnit/TileDB-sumstat/blob/main/CONTRIBUTING.md) on how to contribute to this repo.

#### Additional link to help you get to know TileDB more

- [TileDB](https://tiledb.com/)

---
### Installation

To install this program run the following:

#### Create a virtual environment (recommended)
conda create --name tdbsumstat --file conda-linux-64.lock

conda activate tdbsumstat

#### Install the package
make install

#### Verify installation
tdbsumstat --help

#### Run with singularity

Alternatively you can use the singularity container that already pack tdbsumstat. you can find the container here: ghcr.io/bruno-ariano/tdbsumstat:1.0


#### Run with nextflow

If you are running the pipeline with nextflow change the path where the conda env is stored in conf/base.config

#### To check how to use the Nextflow pipeline for extracting and ingesting data please refer here [HERE](https://github.com/HTGenomeAnalysisUnit/TileDB-sumstat/blob/nextflow_branch/docs/README.md)



