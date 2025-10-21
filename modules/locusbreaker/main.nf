#!/usr/bin/env nextflow

process EXPORT_LOCUSBREAKER {
  label "process_multi"
  conda '/ssu/gassu/conda_envs/scqtl'

  publishDir "${params.outdir}/gwas_and_loci_tables", mode: params.publish_dir_mode


// Define input
  input:
  tuple  val(batch_index), path(traits_list_table)

// Define output
  output:
    path("*_interval.csv"), emit:locus_breaker_tdb_intervals, optional: true
    path("*_segment.csv"), emit:locus_breaker_tdb_segments, optional: true

// Define the shell script to execute
  script:
    """
    scqtl export \
      --table-lb ${traits_list_table} \
      --tiledb-path ${params.tiledb_path} \
      --out ${params.out} \
      --maf-lb ${params.maf_lb} \
      --type-sumstat ${params.type_sumstat} \
      --hole-lb ${params.hole_lb} \
      --locus-max-size-lb ${params.locus_max_size_lb} \
      --locusbreaker \
      --batch-name-lb ${batch_index}
    """
}
