params {
            ingest = true
            file_path_ingestion = "/nfs/users/nfs_b/ba13/TileDB-sumstat/example_data/trait_list.csv"
            mapping_file = "/nfs/users/nfs_b/ba13/TileDB-sumstat/example_data/mapping_file_test.csv"
            ingestion = true
            type_sumstat = "qtl"
            mac = 10
            qc = false
            ingestion_chunk_files = 1
            tiledb_name = "TileDB_test"
            outdir = "/nfs/users/nfs_b/ba13/TileDB-sumstat/test"  // Add this missing parameter
            publish_dir_mode = 'copy'  // Add this missing parameter
        }
