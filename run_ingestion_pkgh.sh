#!/bin/bash
#BSUB -n 2
#BSUB -M 4G
#BSUB -q normal
#BSUB -W 12:00 # time in HH:MM - don't put seconds!
#BSUB -G team151
#BSUB -R "select[mem>4G] rusage[mem=4G] span[hosts=1]"
#BSUB -o tiledb_ingestion.pkgh.log  # Log file for each job
#BSUB -e tiledb_ingestion.pkgh.err  # Error file for each job
set -eo pipefail

#Add the following line if you need to activate conda envs in your script
module load HGI/common/conda
source activate /software/cardinal_analysis/ht/conda_envs/tdbsumstat
module load HGI/common/nextflow/25.04.6
nextflow run ../TileDB-sumstat/main.nf -c ingestion_config_pkgh.nf -profile sanger,conda -work-dir workdir_tiledb -resume
