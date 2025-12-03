process RECOMPUTE_META{
    label 'process_high'
    publishDir "${params.outdir}/gwas_and_loci_tables/", mode: params.publish_dir_mode


        // Define input
    input:
    tuple  val(batch_index), path(traits_list_table)
  
    // Define output
    output:
        path("*.csv"), emit:ltbd_traits, optional: true
    
    // Define the shell script to execute
    script:
        """
        tdbsumstat \
        export \
        --recompute-meta \
        --mac ${params.mac} \
        --trait-list ${traits_list_table} \
        --uri-path ${params.uri_path} \
        --batch-name ${batch_index} \
        --type-sumstat ${params.type_sumstat}
        """
}
