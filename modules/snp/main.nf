process EXPORT_SNP {
    label "process_high"
    publishDir "${params.outdir}/snp_table", mode: params.publish_dir_mode

    // Define input
    input:
    tuple val(chr), val(rows)

    // Define output
    output:
    path("${params.out}_*.csv"), emit: snp_tdb_positions, optional: true

    // Define the shell script to execute
    script:
    // Create the complete CSV content as a Groovy string
    def csvContent = "CHR,POS,TRAIT\n" + rows.collect { "${it.CHR},${it.POS},${it.TRAIT}" }.join('\n')
    
    """
    # Write the complete CSV file
    cat << 'EOF' > batch_${chr}.csv
${csvContent}
EOF
    
    # Run the tdbsumstat command
    tdbsumstat export \
      --snp batch_${chr}.csv \
      --attr ${params.attrs} \
      --tiledb-path ${params.tiledb_path} \
      --out ${params.out} \
      --type-sumstat ${params.type_sumstat}
    """
}