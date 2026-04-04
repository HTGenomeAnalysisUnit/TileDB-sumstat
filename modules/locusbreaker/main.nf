#!/usr/bin/env nextflow

process EXPORT_LOCUSBREAKER {
    label "process_high"
    publishDir "${params.outdir}/gwas_and_loci_tables", mode: params.publish_dir_mode

    // Define input
    input:
    tuple val(batch_index), path(traits_list_table)

    // Define output
    output:
    path("*_interval.csv"), emit: locus_breaker_tdb_intervals, optional: true
    path("*_segment.csv"), emit: locus_breaker_tdb_segments, optional: true

    // Define the shell script to execute
    script:
    """
    tdbsumstat export \
      --table-lb ${traits_list_table} \
      --uri-path ${params.uri_path} \
      --out ${params.out} \
      --maf-lb ${params.maf_lb} \
      --type-sumstat ${params.type_sumstat} \
      --hole-lb ${params.hole_lb} \
      --locus-max-size-lb ${params.locus_max_size_lb} \
      --locusbreaker \
      --batch-name ${batch_index}
    """

    stub:
    """
    touch stub_${batch_index}_interval.csv
    touch stub_${batch_index}_segment.csv
    """
}
