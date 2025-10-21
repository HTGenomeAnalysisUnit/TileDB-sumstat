import click
import cloup
import os
from scqtl.cli.ingestion import ingest
from scqtl.cli.export import export


@cloup.group(name="scqtl", help="Single cell analysis with TileDB", no_args_is_help=True)
@click.pass_context
def cli_init(ctx):
    """Initialize the scqtl CLI context."""
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