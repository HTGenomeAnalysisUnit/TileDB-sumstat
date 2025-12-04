import click
import cloup
import os
from tdbsumstat.cli.ingestion import ingest
from tdbsumstat.cli.export import export


@cloup.group(name="tdbsumstat", help="Genomics summary statistics analysis with TileDB", no_args_is_help=True)
@click.pass_context
def cli_init(ctx):
    """Initialize the tdbsumstat CLI context."""
    ctx.obj = {}
    return ctx.obj


def main():
    # Register subcommands
    cli_init.add_command(ingest)
    cli_init.add_command(export)

    # Invoke the CLI
    cli_init(obj={})


if __name__ == "__main__":
    main()