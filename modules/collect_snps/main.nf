process COLLECT_SNPS {

lebel "process_low"
publishDir "${params.outdir}/snp_table", mode: params.publish_dir_mode

//Define input
input:
path(snp_file)

output:
path("snps_append.csv")

script:

"""
head -n 1 ${snp_file}[0] > all_snps.csv
for f in ${csvs[@]}; do
        tail -n +2 "$f" >> all_snps.csv
    done
"""

}
