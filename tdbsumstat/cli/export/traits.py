"""Trait-based bulk export from TileDB."""
import pandas as pd
import tiledb


def export_by_traits(
    uri_path: str,
    trait_list: str,
    attr: str,
    type_sumstat: str,
    out: str,
    batch_name: str,
) -> None:
    """Export all summary statistics for a list of traits/cells-genes to CSV.

    The function streams data in chunks to handle large datasets efficiently.

    Parameters
    ----------
    uri_path : str
        Path to the TileDB array.
    trait_list : str
        Path to a CSV file with a TRAIT column (and optionally CELL/GENE columns for QTL).
    attr : str
        Comma-separated list of attributes to export.
    type_sumstat : str
        Type of summary statistics: "gwas" or "qtl".
    out : str
        Output file prefix.
    batch_name : str
        Batch identifier appended to the output file name.
    """
    trait_list_pd = pd.read_csv(trait_list)

    with tiledb.open(uri_path, mode="r") as array:
        if type_sumstat == "gwas":
            trait_list_np = trait_list_pd["TRAIT"].to_list()
            tiledb_iterator = array.query(
                return_incomplete=True,
                attrs=attr.split(","),
            ).df[:, trait_list_np, :]
        else:
            trait_list_pd[["cell", "gene"]] = trait_list_pd["TRAIT"].str.split(":", expand=True)
            cells = trait_list_pd["cell"].to_list()
            gene = trait_list_pd["gene"].to_list()
            tiledb_iterator = array.query(
                return_incomplete=True,
                attrs=attr.split(","),
            ).df[:, cells, gene, :]

        for chunk in tiledb_iterator:
            chunk.to_csv(f"{out}_{batch_name}.csv", mode="a", index=False, header=False)

    print(f"Saved filtered summary statistics in {out}")
