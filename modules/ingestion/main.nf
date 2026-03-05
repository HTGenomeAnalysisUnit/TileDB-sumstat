#!/usr/bin/env nextflow

process INGEST_DATA {
    label "process_high"
    //publishDir "${params.outdir}/TileDB/", mode: 'link'

    // Define input
    input:
    path(tiledb)
    each path(list_files)
    path(mapping_file)
    path(dummy_file)

    // Define output
    output:
    path("TileDB_${params.tiledb_name}_metadata_parts_${list_files.name}"), emit: metadata_parts
    path("ingestion_complete_${list_files.name}"), emit: ingestion_done
    path("${tiledb}"), emit: tiledb_updated

    // Define the shell script to execute
    script:
    def qc = params.qc ? "--qc" : ""
    def pvar_file = params.pvar_file ? "--pvar-file ${params.pvar_file}" : ""
    def permuted = params.permuted ? "--permuted" : ""
    """
    tdbsumstat ingest \
    --uri-path TileDB_${params.tiledb_name}\
    --file-path ${list_files} \
    --mapping-file ${mapping_file} \
    --type-sumstat ${params.type_sumstat} \
    --maf ${params.maf} \
    --mac ${params.mac} ${qc} ${pvar_file} ${permuted}

    # Rename the metadata parts directory to include the list_files name for uniqueness
    if [ -d "TileDB_${params.tiledb_name}_metadata_parts" ]; then
        mv "TileDB_${params.tiledb_name}_metadata_parts" "TileDB_${params.tiledb_name}_metadata_parts_${list_files.name}"
    else
        # Create empty directory if none was created
        mkdir -p "TileDB_${params.tiledb_name}_metadata_parts_${list_files.name}"
    fi

    # Create a dummy file to signal completion
    touch ingestion_complete_${list_files.name}
    """

    stub:
    """
    mkdir -p TileDB_${params.tiledb_name}_metadata_parts_${list_files.name}
    touch ingestion_complete_${list_files.name}
    """
}
