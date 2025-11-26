#!/usr/bin/env nextflow

process CREATE_TILEDB {
    label "process_high"
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
    tdbsumstat ingest --uri-path TileDB_${params.uri_path} \
    --create-tiledb \
    --mapping-file ${mapping_file} \
    --type-sumstat ${params.type_sumstat}

    touch dummy_file
    """
}