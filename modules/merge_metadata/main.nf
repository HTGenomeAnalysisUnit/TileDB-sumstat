#!/usr/bin/env nextflow

process MERGE_METADATA {
    label "process_high"
    publishDir "${params.outdir}/TileDB/", mode: 'copy', pattern: "TileDB_${params.tiledb_name}"
    publishDir "${params.outdir}/TileDB/", mode: 'copy', pattern: "TileDB_${params.tiledb_name}_metadata.csv"

    input:
    path(tiledb)
    path(mapping_file)
    path(metadata_parts)
    path(ingestion_signals)

    output:
    path("${tiledb}"), emit: tiledb_final
    path("TileDB_${params.tiledb_name}_metadata.csv"), emit: metadata
    path("tiledb_with_metadata"), emit: completion_signal

    script:
    """
    # Create the main metadata_parts directory
    mkdir -p TileDB_${params.tiledb_name}_metadata_parts

    # Copy all JSON files from all input directories
    for parts_dir in ${metadata_parts}; do
      if [ -d "\$parts_dir" ]; then

        # Extract the trailing number before .csv from the input folder
        base=\$(basename "\$parts_dir")
        base="\${base%.csv}"
        suffix="\${base##*.}"

        for json in "\$parts_dir"/*.json; do
          [ -e "\$json" ] || continue

          # Add that trailing number to .json file name
          base=\$(basename "\$json")
          name="\${base%.json}"

          cp -v "\$json" "TileDB_${params.tiledb_name}_metadata_parts/\${name}_\${suffix}.json"
        done
      fi
    done

    echo "Total JSON files collected:"
    ls TileDB_${params.tiledb_name}_metadata_parts/*.json 2>/dev/null | wc -l || echo "0"

    tdbsumstat ingest \
    --uri-path ${tiledb} \
    --mapping-file ${mapping_file} \
    --type-sumstat ${params.type_sumstat} \
    --only-meta

    # Signal that metadata merge is complete
    touch tiledb_with_metadata
    """

    stub:
    """
    touch TileDB_${params.tiledb_name}_metadata.csv
    touch tiledb_with_metadata
    """
}
