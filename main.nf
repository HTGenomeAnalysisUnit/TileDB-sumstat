include { INGEST_DATA } from "./modules/ingestion"
include { CREATE_TILEDB } from "./modules/create_tiledb"
include { EXPORT_LOCUSBREAKER } from "./modules/locusbreaker"
include { EXPORT_SNP } from "./modules/snp"
include { EXPORT_REGIONS } from "./modules/regions"
include { MERGE_METADATA } from "./modules/merge_metadata"
include { TRAITS } from "./modules/traits"
include { RECOMPUTE_META } from "./modules/recompute_meta"

workflow {
    if (params.ingestion){
        Channel.fromPath(params.file_path_ingestion, checkIfExists:true)
        .splitText(by: params.ingestion_chunk_files, keepHeader: true, file: true)
        .set { list_files }
        
        mapping_file = Channel.fromPath(params.mapping_file, checkIfExists:true)
        
        // Check if TileDB already exists
        def tiledb_path = "${params.outdir}/TileDB/${params.tiledb_name}"
        def tiledb_exists = file(tiledb_path).exists()
        
        if (!tiledb_exists) {
            // Create TileDB if it doesn't exist
            create_tiledb = CREATE_TILEDB(mapping_file, Channel.of('dummy'))
            tiledb_storage = create_tiledb.tiledb_storage
            dummy_file = create_tiledb.dummy_file
        } else {
            // Use existing TileDB
            tiledb_storage = Channel.fromPath(tiledb_path)
            dummy_file = Channel.fromPath("${tiledb_path}/.dummy").ifEmpty { 
                Channel.of(file("${tiledb_path}/.dummy").tap { it.text = "dummy" }) 
            }
        }
        
        // Pass the TileDB array through all ingestion steps
        ingestion_results = INGEST_DATA(tiledb_storage, list_files, mapping_file, dummy_file)
        
        // Rest of your workflow...
        all_metadata_parts = ingestion_results.metadata_parts.collect()
        all_ingestion_done = ingestion_results.ingestion_done.collect()
        updated_tiledb = ingestion_results.tiledb_updated.first()
        merged_metadata = MERGE_METADATA(updated_tiledb, mapping_file, all_metadata_parts, all_ingestion_done) 
    }
    if (params.export){
        if (params.snp) {
        Channel
        .fromPath(params.snp, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            tuple(row.CHR, row)
        }
        .groupTuple()
        .set { tiledb_metadata_batches }
        snp_results = EXPORT_SNP(tiledb_metadata_batches)
	snp_results.snp_tdb_positions
        	.transpose()  // Convert tuple(chr, [file1, file2]) to [tuple(chr, file1), tuple(chr, file2)]
        	.collectFile(keepHeader:true) { chr, file ->
            	["${params.out}_chr_${chr}_concatenated.csv", file.text]
        	}
        	.set { concatenated_files }
    
    // Optionally publish the concatenated files
    concatenated_files.subscribe { file ->
        file.copyTo("${params.outdir}/snp_table_concatenated/${file.name}")
    }
        }
        }
        if (params.locusbreaker){
            Channel.fromPath(params.table_lb, checkIfExists:true)
                .splitText(by: params.tiledb_batch_size, keepHeader: true, file: true)
                .map { batch_file -> 
                def batch_index = (batch_file.name =~ /\.(\d+)\.csv$/)[0][1]
                tuple(batch_index, batch_file)
                }.set { tiledb_metadata_batches }
                
            EXPORT_LOCUSBREAKER(tiledb_metadata_batches)
         }
        if (params.export_traits){
            Channel.fromPath(params.list_traits, checkIfExists:true)
                .splitText(by: params.tiledb_batch_size, keepHeader: true, file: true)
                .map { batch_file -> 
                def batch_index = (batch_file.name =~ /\.(\d+)\.csv$/)[0][1]
                tuple(batch_index, batch_file)
                }.set { tiledb_metadata_batches }

            TRAITS(tiledb_metadata_batches)
        
    }
        if (params.recompute_meta){
            Channel.fromPath(params.list_traits, checkIfExists:true)
                .splitText(by: params.tiledb_batch_size, keepHeader: true, file: true)
                .map { batch_file -> 
                def batch_index = (batch_file.name =~ /\.(\d+)\.csv$/)[0][1]
                tuple(batch_index, batch_file)
                }.set { tiledb_metadata_batches }

            RECOMPUTE_META(tiledb_metadata_batches)
        
    }
}
