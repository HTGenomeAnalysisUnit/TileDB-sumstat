import pandas as pd
import numpy as np
import click

@click.command()
@click.option("--num_snps", default=10000, help="Total number of SNPs to create")
@click.option("--num_snps_gene", default=2000, help="How many SNPs each gene should have")
@click.option("--out_csv", default="dummy_out", help="Where to send the output")

def create_dummy_data(num_snps, num_snps_gene, out_csv):
    # Initialize lists for each column

    gene_ids = []
    variant_ids = []
    start_distances = []
    afs = []
    ma_samples = []
    ma_counts = []
    pval_nominals = []
    slopes = []
    slope_ses = []
    positions = []
    a1s = []
    a2s = []
    
    # Keep track of generated variant IDs to avoid duplicates
    generated_variant_ids = set()
    num_genes = 2

    #num_genes = round(num_snps / num_snps_gene)
    for gene_idx in range(num_genes):
        # Generate unique phenotype_id for this gene
        phenotype_id = f"ENSG00000{num_snps + gene_idx}"
        for _ in range(num_snps_gene):
            # Ensure unique variant_id
            while True:
                pos = np.random.randint(1, num_snps)
                a1 = np.random.choice(['A', 'T', 'C', 'G'])
                a2 = np.random.choice(['A', 'T', 'C', 'G'])
                while a1==a2:
                    a1 = np.random.choice(['A', 'T', 'C', 'G'])
                    a2 = np.random.choice(['A', 'T', 'C', 'G'])
                
                variant_id = f"chr20:{pos}:{a1}:{a2}"
                variant_reverse_id = f"chr20:{pos}:{a2}:{a1}"

                if variant_id not in generated_variant_ids and variant_reverse_id not in generated_variant_ids:
                    generated_variant_ids.add(variant_id)
                    generated_variant_ids.add(variant_reverse_id)
                    break
                
            
            # Generate other columns with dummy data
            start_distance = np.random.randint(-500000, 500000)
            a1s.append(a1)
            a2s.append(a2)
            af = np.random.uniform(0, 1)
            ma_sample = np.random.randint(50, 200)
            ma_count = np.random.randint(50, 300)
            pval_nominal = np.random.uniform(0, 1)
            slope = np.random.uniform(-1, 1)
            slope_se = np.random.uniform(0.001, 1)

            # Append data to lists
            gene_ids.append(phenotype_id)
            variant_ids.append(variant_id)
            start_distances.append(start_distance)
            afs.append(af)
            positions.append(pos)
            ma_samples.append(ma_sample)
            ma_counts.append(ma_count)
            pval_nominals.append(pval_nominal)
            slopes.append(slope)
            slope_ses.append(slope_se)

        # Create a DataFrame
        data = {
        "Chr":20,
        "Gene":gene_ids,
        "cell.type": "Tgd",
        "pos": positions,
        "a0":a1,
        "a1":a2,
        "SNP": variant_ids,
        "START": start_distances,
        "EAF": afs,
        "p": pval_nominals,
        "beta": slopes,
        "se": slope_ses
        }

        df = pd.DataFrame(data)

        # Save to a CSV file
        df.to_csv(f"{out_csv}_{phenotype_id}.tsv.gz", index=False, sep="\t", compression="gzip")

        print(f"Dummy data file created: {out_csv}")

if __name__ == '__main__':
    create_dummy_data()
