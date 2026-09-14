process TRAITS{
    label 'process_high'
    publishDir "${params.outdir}/gwas_and_loci_tables/", mode: params.publish_dir_mode


        // Define input
    input:
    tuple  val(batch_index), path(traits_list_table)
  
    // Define output
    output:
        path("*_${batch_index}.csv"), emit:ltbd_traits, optional: true
    
    // Define the shell script to execute
    script:
        // params.out is a file name prefix, not a directory: only its last component is
        // kept so that results are written in the task dir and collected by publishDir
        def out_name = params.out ? params.out.toString().replaceAll('.*/', '') : ''
        def out_prefix = out_name in ['', '.'] ? 'traits' : out_name
        """
        tdbsumstat \
        export \
        --trait-list ${traits_list_table} \
        --uri-path ${params.uri_path} \
        --attr ${params.attrs} \
        --out ${out_prefix} \
        --batch-name ${batch_index} \
        --type-sumstat ${params.type_sumstat}
        """
}
