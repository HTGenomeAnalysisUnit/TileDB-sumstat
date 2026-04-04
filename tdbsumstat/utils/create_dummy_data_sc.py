"""Generate synthetic per-gene QTL tables for local tests (alleles and SNP column aligned)."""
from __future__ import annotations

import argparse
import csv
import gzip
import random
from typing import Optional

# Base genomic position offset per gene so (CHR, CELL, GENE, POS) never collides across genes
# when ingested into the same TileDB array.
_POS_OFFSET_PER_GENE = 100_000

_ALLELES = ("A", "T", "C", "G")

_HEADER = [
    "Chr",
    "Gene",
    "cell.type",
    "pos",
    "a0",
    "a1",
    "SNP",
    "START",
    "EAF",
    "p",
    "beta",
    "se",
]


def _build_variant_key(chrom: int, pos: int, al1: str, al2: str) -> tuple[str, str]:
    """Return ordered pair of variant id strings for a locus (forward and reverse).

    Used to deduplicate the same physical SNP when alleles are swapped.

    Example:
        >>> _build_variant_key(20, 100, "A", "G")
        ('chr20:100:A:G', 'chr20:100:G:A')
    """
    fwd = f"chr{chrom}:{pos}:{al1}:{al2}"
    rev = f"chr{chrom}:{pos}:{al2}:{al1}"
    return fwd, rev


def _random_allele_pair(rng: random.Random) -> tuple[str, str]:
    """Draw two distinct nucleotides uniformly from A/T/C/G.

    Example:
        >>> _random_allele_pair(random.Random(0))
        ('T', 'C')
    """
    while True:
        a1 = rng.choice(_ALLELES)
        a2 = rng.choice(_ALLELES)
        if a1 != a2:
            return a1, a2


def _generate_one_gene_rows(
    rng: random.Random,
    num_snps: int,
    num_snps_gene: int,
    gene_idx: int,
    generated_variant_ids: set[str],
) -> tuple[str, list[dict]]:
    """Build table rows for a single gene; returns (phenotype_id, rows).

    Example:
        >>> r = random.Random(0)
        >>> g, rows = _generate_one_gene_rows(r, 100, 3, 0, set())
        >>> len(rows)
        3
    """
    phenotype_id = f"ENSG00000{num_snps + gene_idx}"
    pos_offset = gene_idx * _POS_OFFSET_PER_GENE
    rows: list[dict] = []

    for _ in range(num_snps_gene):
        while True:
            pos_local = rng.randrange(1, num_snps)
            pos = pos_local + pos_offset
            a1, a2 = _random_allele_pair(rng)

            variant_id, variant_reverse_id = _build_variant_key(20, pos, a1, a2)

            if variant_id not in generated_variant_ids and variant_reverse_id not in generated_variant_ids:
                generated_variant_ids.add(variant_id)
                generated_variant_ids.add(variant_reverse_id)
                break

        rows.append(
            {
                "Chr": 20,
                "Gene": phenotype_id,
                "cell.type": "Tgd",
                "pos": pos,
                "a0": a1,
                "a1": a2,
                "SNP": variant_id,
                "START": rng.randrange(-500_000, 500_000),
                "EAF": rng.random(),
                "p": rng.random(),
                "beta": rng.uniform(-1, 1),
                "se": rng.uniform(0.001, 1),
            }
        )

    return phenotype_id, rows


def create_dummy_data(
    num_snps: int,
    num_snps_gene: int,
    out_csv: str,
    seed: Optional[int],
    also_tsv: bool,
) -> None:
    """Write two gzipped TSVs with disjoint positions per gene, random A/T/C/G alleles, SNP = chr:pos:a0:a1.

    Example:
        ``create_dummy_data(10000, 2000, "dummy_out", seed=42, also_tsv=True)``
        writes ``dummy_out_ENSG0000010000.tsv.gz`` (and optionally ``.tsv``) with varied alleles per row.
    """
    rng = random.Random(seed)
    num_genes = 2
    generated_variant_ids: set[str] = set()

    for gene_idx in range(num_genes):
        phenotype_id, rows = _generate_one_gene_rows(
            rng, num_snps, num_snps_gene, gene_idx, generated_variant_ids
        )

        gz_path = f"{out_csv}_{phenotype_id}.tsv.gz"
        with gzip.open(gz_path, "wt", encoding="utf-8", newline="") as gz_f:
            w = csv.DictWriter(gz_f, fieldnames=_HEADER, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        print(f"Dummy data file created: {gz_path}")

        if also_tsv:
            tsv_path = f"{out_csv}_{phenotype_id}.tsv"
            with open(tsv_path, "w", encoding="utf-8", newline="") as tsv_f:
                w = csv.DictWriter(tsv_f, fieldnames=_HEADER, delimiter="\t", lineterminator="\n")
                w.writeheader()
                w.writerows(rows)
            print(f"Uncompressed copy: {tsv_path}")


def _parse_args() -> argparse.Namespace:
    """Parse CLI arguments for the dummy-data generator."""
    p = argparse.ArgumentParser(description="Generate two dummy QTL TSVs (random A/T/C/G alleles per row).")
    p.add_argument("--num_snps", type=int, default=10000, help="Upper bound for random POS (exclusive)")
    p.add_argument("--num_snps_gene", type=int, default=2000, help="SNPs per gene file")
    p.add_argument("--out_csv", type=str, default="dummy_out", help="Output prefix")
    p.add_argument("--seed", type=int, default=None, help="RNG seed (reproducible runs)")
    p.add_argument("--also-tsv", action="store_true", help="Also write uncompressed .tsv files")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    create_dummy_data(
        num_snps=args.num_snps,
        num_snps_gene=args.num_snps_gene,
        out_csv=args.out_csv,
        seed=args.seed,
        also_tsv=args.also_tsv,
    )
