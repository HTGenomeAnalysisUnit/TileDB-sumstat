#!/usr/bin/env nextflow

process LOCUS_BREAKER_TILEDB {
  label "process_multi"
  // conda '/ssu/gassu/conda_envs/scqtl'
  conda '/software/cardinal_analysis/ht/conda_envs/scqtl'

  publishDir "${params.outdir}/results/gwas_and_loci_tables/", mode: params.publish_dir_mode


// Define input
  input:
  tuple  path(traits_list_table_ingest)

// Define output
  output:
    path("tiledb_uri"), emit:tiledb_storage, optional: true
    path("metadata.csv"), emit:metadata

// Define the shell script to execute
  script:
    """
    scqtl --workers ${params.workers}  --memory_w ${params.memory_w} \
    ingest
    --uri-path ${params.tiledb_uri} \ 
    --file-path ${traits_list_table_ingest} \
    --mapping-file ${mapping_file} \
    --type-sumstat ${params.type_sumstat} \
    --chunk-files ${meta}\
    """
}
