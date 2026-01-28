import os
import gwaslab as gl
import polars as pl
import pandas as pd

def _qc_sumstat(
                uri: str = None,
                type_sumstat: str = None,
                type_trait:str = None,
                file_path:str = None):
        directory = uri + "_logs"
        filename = os.path.basename(file_path)   # "test.csv.gz"
        # Remove all extensions
        file_name = filename.split('.')[0]       # "test"
    
        if not os.path.isdir(directory):
            os.mkdir(directory)
        sumstat_preqc = chunk_pl.to_pandas()
        if type_sumstat == "gwas":
            if type_trait== "quant":
                sumstat_gl =gl.Sumstats(sumstat_preqc,
                    snpid="SNPID",
                    chrom="CHR",
                    pos="POS",
                    eaf="EAF",
                    beta="BETA",
                    se="SE",
                    p="P",
                    n="N",
                    ea = "EA",
                    nea = "NEA",
                    other = ["TRAIT","RSID"])
            else:
                sumstat_gl =gl.Sumstats(sumstat_preqc,
                    snpid="SNPID",
                    chrom="CHR",
                    pos="POS",
                    eaf="EAF",
                    beta="BETA",
                    se="SE",
                    p="P",
                    n="N",
                    ncase = "N_CASES",
                    ncontrol = "N_CONTROLS",
                    ea = "EA",
                    nea = "NEA",
                    other = ["TRAIT","RSID"])

        else:
            sumstat_gl =gl.Sumstats(sumstat_preqc,
                 snpid="SNPID",
                 chrom="CHR",
                 pos="POS",
                 eaf="EAF",
                 beta="BETA",
                 se="SE",
                 p="P",
                 n="N",
                 other = ["CELL","GENE","RSID","DIST","PHENO_VAR"])
        #sumstat_gl.fix_id()
        sumstat_gl.fix_chr(remove=True)
        sumstat_gl.fix_pos(remove=True)
        sumstat_gl.fix_allele(remove=True)
        sumstat_gl.check_sanity()
        sumstat_gl.check_data_consistency()
        #sumstat_gl.remove_dup(mode="m")
        #sumstat_gl.basic_check(n_cores = 4, remove=True, remove_dup=True)

        sumstat_gl.log.save(directory + "/" + file_name)
        chunk_pl = pl.from_pandas(sumstat_gl.data)