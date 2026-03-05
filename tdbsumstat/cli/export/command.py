"""CLI command definition for the ``export`` subcommand.

This module wires CLI options to the individual export handler functions
defined in the sibling modules (snp, regions, locusbreaker, metadata, traits).
"""
import click
import cloup

from tdbsumstat.cli.export.helpers import open_tiledb_and_load_metadata
from tdbsumstat.cli.export.locusbreaker import export_with_locusbreaker
from tdbsumstat.cli.export.metadata import export_metadata, recompute_metadata
from tdbsumstat.cli.export.regions import export_by_regions
from tdbsumstat.cli.export.snp import export_by_snp
from tdbsumstat.cli.export.traits import export_by_traits

help_doc = """
Query TileDB database and export data.
"""


@cloup.command("export", no_args_is_help=True, help=help_doc)
@cloup.option_group(
    "Options for querying specific chromosomes, cells, genes or positions in the TileDB",
    cloup.option("--uri-path", default=None, type=str, help="path of TileDB"),
    cloup.option("--table-regions", default=None, type=str, help="Regions to interrogate from a table"),
    cloup.option("--trait-list", default=None, type=str, help="List of entire traits to filter"),
    cloup.option("--attr", default="P,SNPID,EAF,BETA,SE", type=str, help="Attributes to output"),
    cloup.option("--export-meta", is_flag=True, default=False, type=str, help="Get metadata from TileDB"),
    cloup.option("--mac", default=0, type=int, help="Filter for MAC when recomputing metadata"),
    cloup.option(
        "--recompute-meta",
        is_flag=True,
        default=False,
        type=str,
        help="Recompute metadata after applying filters (Does not modify data within the TileDB)",
    ),
    cloup.option(
        "--snp",
        default=None,
        type=str,
        help="List of SNPs to interrogate taken from a txt file. Please check README for details on the format of this file",
    ),
    cloup.option("--batch-name", default=None, type=str, help="Name of the batch"),
)
@cloup.option_group(
    "Options for Locusbreaker",
    cloup.option("--locusbreaker", is_flag=True, type=bool, default=False, help="Option to run locusbreaker"),
    cloup.option(
        "--hole-lb",
        default=250000,
        type=int,
        help="Minimum base-pair distance between SNPs in different loci (default: 250000)",
    ),
    cloup.option(
        "--maf-lb",
        default=0.001,
        type=float,
        help="The MAF to filter the TILEDB before locusbreaker is run",
    ),
    cloup.option(
        "--locus-max-size-lb",
        default=3000000,
        type=float,
        help="The maximum size allowed for the locus. Default: 1Mb",
    ),
    cloup.option(
        "--cis-trans-lb",
        default="cis",
        type=str,
        help="If locusbreaker run on cis or trans QLTs",
    ),
    cloup.option("--table-lb", default=None, type=str, help="Path of the table to provide"),
    cloup.option("--type-sumstat", default=None, type=str, help="Type of summary data"),
)
@cloup.option_group(
    "Options for output",
    cloup.option(
        "--out",
        default="out",
        type=str,
        help="Output path with file name where results will be stored",
    ),
)
@click.pass_context
def export(
    ctx,
    uri_path: str,
    type_sumstat: str,
    table_regions: str,
    trait_list: str,
    mac: int,
    attr: str,
    snp: str,
    export_meta: bool,
    recompute_meta: bool,
    locusbreaker: bool,
    maf_lb: float,
    cis_trans_lb: str,
    table_lb: str,
    hole_lb: int,
    out: str,
    locus_max_size_lb: int,
    batch_name: str,
):
    """Route the export command to the appropriate handler based on CLI flags."""
    if snp:
        tiledb_export, _df_meta = open_tiledb_and_load_metadata(uri_path, type_sumstat)
        export_by_snp(tiledb_export, snp, attr, type_sumstat, out)
        tiledb_export.close()

    elif table_regions:
        tiledb_export, _df_meta = open_tiledb_and_load_metadata(uri_path, type_sumstat)
        export_by_regions(tiledb_export, table_regions, attr, type_sumstat, out)
        tiledb_export.close()

    elif locusbreaker:
        _tiledb_export, df_meta = open_tiledb_and_load_metadata(uri_path, type_sumstat)
        _tiledb_export.close()
        export_with_locusbreaker(
            uri_path=uri_path,
            df_meta=df_meta,
            table_lb=table_lb,
            maf_lb=maf_lb,
            hole_lb=hole_lb,
            locus_max_size_lb=locus_max_size_lb,
            cis_trans_lb=cis_trans_lb,
            type_sumstat=type_sumstat,
            out=out,
            batch_name=batch_name,
        )

    elif export_meta:
        export_metadata(uri_path, type_sumstat, out)

    elif recompute_meta:
        tiledb_export, df_meta = open_tiledb_and_load_metadata(uri_path, type_sumstat)
        recompute_metadata(tiledb_export, df_meta, trait_list, type_sumstat, mac, out, batch_name)
        tiledb_export.close()

    else:
        export_by_traits(uri_path, trait_list, attr, type_sumstat, out, batch_name)
