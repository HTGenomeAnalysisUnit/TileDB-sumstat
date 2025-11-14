#!/usr/bin/env nextflow

process EXPORT_SNP {

  publishDir "${params.outdir}/snp_table", mode: params.publish_dir_mode


// Define input
  input:
  tuple  val(batch_index), path(snps_list_table)

// Define output
  output:
    path("*_batch_*.csv"), emit:snp_tdb_positions, optional: true

// Define the shell script to execute
  script:
    """
    tdbsumstat export \
      --snp ${snps_list_table} \
      --attr ${params.attrs} \
      --tiledb-path ${params.tiledb_path} \
      --out ${params.out} \
      --type-sumstat ${params.type_sumstat} \
      --batch-name ${batch_index}
    """
}
