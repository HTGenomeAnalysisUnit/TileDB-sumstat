process EXPORT_SNP {
    label "process_high"

    // Define input
    input:
    tuple val(chr), val(rows)

    // Define output
    output:
    tuple val(chr), path("${params.out}*.csv"), emit: snp_tdb_positions, optional: true

    // Define the shell script to execute
    script:
    // Create the complete CSV content as a Groovy string
    // Get column names and create header
    def columnNames = rows[0].keySet().toList().sort()
    def header = columnNames.join(',')
    
    """
    # Create header
    echo "${header}" > snp_batch_${chr}.csv
    # Use collectFile-like approach with shell commands
    ${rows.collect { row ->
        def line = columnNames.collect { col -> row[col] ?: '' }.join(',')
        "echo '${line}' >> snp_batch_${chr}.csv"
    }.join('\n    ')}    
    # Run the tdbsumstat command
    tdbsumstat export \
      --snp snp_batch_${chr}.csv \
      --attr ${params.attrs} \
      --tiledb-path ${params.tiledb_path} \
      --out ${params.out} \
      --type-sumstat ${params.type_sumstat}
    """
}
