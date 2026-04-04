#!/usr/bin/env nextflow

process EXPORT_REGIONS {
    label "process_high"
    publishDir "${params.outdir}/results/regions_export/", mode: params.publish_dir_mode

    // Define input
    input:
    tuple val(batch_index), path(regions_table)

    // Define output
    output:
    path("*.csv"), emit: regions_output, optional: true

    // Define the shell script to execute
    script:
    """
    tdbsumstat export \
      --table-regions ${regions_table} \
      --uri-path ${params.uri_path} \
      --type-sumstat ${params.type_sumstat} \
      --out ${params.out}_${batch_index} \
      --attr ${params.attrs}

    # Ensure output has .csv extension
    for f in ${params.out}_${batch_index}*; do
        [ -f "\$f" ] && [[ "\${f}" != *.csv ]] && mv "\$f" "\${f}.csv"
    done
    """

    stub:
    """
    touch stub_regions_${batch_index}.csv
    """
}
