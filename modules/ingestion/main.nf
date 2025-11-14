#!/usr/bin/env nextflow

process INGEST_DATA {
  label "process_multi"
  conda 'tdbsumstat'
  publishDir "${params.outdir}/TileDB/", mode: 'copy'


// Define input
  input:
  path(tiledb)
  each path(list_files)
  path(mapping_file)
  path(dummy_file)

// Define output
  output:
    path("metadata_parts_${list_files.name}"), emit: metadata_parts
    path("ingestion_complete_${list_files.name}"), emit: ingestion_done
    path("${tiledb}"), emit: tiledb_updated

// Define the shell script to execute
  script:
    def qc = params.qc ? "--qc" : ""
    """
    tdbsumstat ingest \
    --uri-path TileDB_${params.tiledb_name}\
    --file-path ${list_files} \
    --mapping-file ${mapping_file} \
    --type-sumstat ${params.type_sumstat} \
    ${qc}

    # Copy all metadata JSON files to a unique directory for output
    mkdir -p metadata_parts_${list_files.name}
    cp -v ${tiledb}_metadata_parts/*.json metadata_parts_${list_files.name}/ || true
  
    # Create a dummy file to signal completion
    touch ingestion_complete_${list_files.name}
    """
}
