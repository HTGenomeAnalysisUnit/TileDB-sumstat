import pandas as pd
import numpy as np
import click

@click.command()
@click.option("--num_snps", default=10000, help="Total number of SNP to create")
@click.option("--num_snps_gene", default=2000, help="How many snps each gene should have")
@click.option("--out_csv", default="dummy_out.tsv.gz", help="Where to send the output")

# Parameters

# Generating the data
def create_dummy_data(num_snps, num_snps_gene, out_csv):
    # Initialize lists for each column
    phenotype_ids = []
    variant_ids = []
    start_distances = []
    afs = []
    ma_samples = []
    ma_counts = []
    pval_nominals = []
    slopes = []
    slope_ses = []
    num_genes = round(num_snps/num_snps_gene)
    for gene_idx in range(num_genes):
        # Generate unique phenotype_id for this gene
        phenotype_id = f"ENSG00000{num_snps + gene_idx}"
        for i in range(num_snps_gene):
            # Generate columns with dummy data
            variant_id = f"chr{np.random.randint(1, 24)}_{np.random.randint(1, num_snps)}_{np.random.choice(['A', 'T', 'C', 'G'])}_{np.random.choice(['A', 'T', 'C', 'G'])}"
            start_distance = np.random.randint(-500000, 500000)
            af = np.random.uniform(0, 1)
            ma_sample = np.random.randint(50, 200)
            ma_count = np.random.randint(50, 300)
            pval_nominal = np.random.uniform(0, 1)
            slope = np.random.uniform(-1, 1)
            slope_se = np.random.uniform(0.001, 1)

            # Append data to lists
            phenotype_ids.append(phenotype_id)
            variant_ids.append(variant_id)
            start_distances.append(start_distance)
            afs.append(af)
            ma_samples.append(ma_sample)
            ma_counts.append(ma_count)
            pval_nominals.append(pval_nominal)
            slopes.append(slope)
            slope_ses.append(slope_se)

    # Create a DataFrame
    data = {
    "phenotype_id": phenotype_ids,
    "variant_id": variant_ids,
    "start_distance": start_distances,
    "af": afs,
    "ma_samples": ma_samples,
    "ma_count": ma_counts,
    "pval_nominal": pval_nominals,
    "slope": slopes,
    "slope_se": slope_ses
    }

    df = pd.DataFrame(data)

    # Save to a CSV file
    df.to_csv(out_csv, index=False, sep = "\t", compression="gzip")

    print(f"Dummy data file created: {out_csv}")

if __name__ == '__main__':
    create_dummy_data()