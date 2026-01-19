library(data.table)
library(tidyr)
library(tiledb)
library(qvalue)
library(pbapply)
library(dplyr)


# Set up cohort
cohort = "UKBB" ###

### Load summary with phenotypic variance - identify also genes expressed in at least 20% of individuals!
pheno_var_files <- Sys.glob(paste0("/lustre/scratch124/humgen/projects_v2/cardinal_analysis/analysed_datasets/eqtls/F3_QTLs_updated_12_2025/TensorQTL_", cohort, "/celltype_2/results_F3_cis/norm_data/dMean__*__all/gene_expression_stats_after_rank_based_INT.tsv"))
pheno_var <- rbindlist(lapply(pheno_var_files, fread)) %>%
	mutate(
		TRAIT = paste0(celltype, ":", gene),
		expressed_at_least_20perc = ((n_total_samples*20)/100) <= n_expressing_samples
	) %>%
	rename(PHENO_VAR=var_expression)

### Load TileDB and metadata
metadata_file <- paste0("/lustre/scratch124/humgen/projects_v2/cardinal_analysis/analysed_datasets/eqtls/F3_QTLs_updated_12_2025/TileDB_ingestion/TileDB_", cohort, "/celltype_2/TileDB/TileDB_tiledb_", tolower(cohort), "_celltype2_f3_16_12_25_metadata.csv")

### Compute qvalues and add correct pheno variance to metadata
metadata <- fread(metadata_file) %>%
	mutate(TRAIT = paste0(CELL, ":", GENE)) %>%
	select(-PHENO_VAR) %>% ### remove original phenotypic variance, set to 1
	left_join(pheno_var %>% select(TRAIT, PHENO_VAR, expressed_at_least_20perc), by="TRAIT") %>%
	filter(expressed_at_least_20perc) %>% ### filter out genes expressed 
	mutate(qvalue = qvalue::qvalue(ACAT)$qvalues)

#Add SIG and LIM values
res_df <- rbindlist(res) %>%
	mutate(
		MIN_P = ifelse(MIN_P==0, 1e-100, MIN_P),
		SIG = MIN_P*2,
		LIM = ifelse(MIN_P>1e-6, MIN_P*10, 1e-5)
	)
# Merge back with original metadata and save
final <- metadata %>%
	filter(qvalue < 0.05) %>%
	select(CHR,TRAIT,SIG,LIM,PHENO_VAR)

fwrite(final, paste0("/lustre/scratch124/humgen/projects_v2/cardinal_analysis/analysis/ba13/finemapping/", cohort, "_cardinal_12_25_finemapping/TileDB_tiledb_", tolower(cohort), "_celltype2_f3_16_12_25_metadata_FINEMAP_INPUT.csv"), sep=",", quote=F, na=NA)

