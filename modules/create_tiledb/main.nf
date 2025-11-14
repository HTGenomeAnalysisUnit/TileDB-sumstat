#!/usr/bin/env nextflow

process CREATE_TILEDB {
    label "process_single"
    conda '/ssu/gassu/conda_envs/tdbsumstat'
    publishDir "${params.outdir}/TileDB/", mode: 'copy'

    input:
    path(mapping_file)
    val dummy
    // Define output
    output:
    path "TileDB_${params.tiledb_name}", emit: tiledb_storage
    path "dummy_file", emit: dummy_file

    // Define the shell script to execute
    script:
    """
    touch dummy_file

    tdbsumstat ingest --uri-path TileDB_${params.tiledb_name} \
    --create-tiledb \
    --mapping-file ${mapping_file} \
    --type-sumstat ${params.type_sumstat}
    """
}