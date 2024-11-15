#Create Dask cluster
from dask_jobqueue import LSFCluster as Cluster
from dask.distributed import Client
import click
import cloup
from cli.ingestion import ingestion
from cli.export import export

@cloup.group(name="main", help="Single cell analysis with TileDB", no_args_is_help=True)
@cloup.option_group(
    "Dask cluster options",
    cloup.option("--workers", default = None, type=int, help = "Number of workers to use in the Dask cluster"),
    cloup.option("--cores_w", default = 1, type=int, help = "Number of core to use per single worker"),
    cloup.option("--memory_w", default = "8", type=str, help = "Number of GB to give as memory for each worker")
)

@click.pass_context
def cli_init(ctx, workers: int, cores_w: int, memory_w: str):

    #Create Dask cluster
    if not workers is None:
        cluster = Cluster(cores = cores_w, memory=f"{memory_w} GB",interface='ib0')
        cluster.scale(workers)
        client = Client(cluster)
        ctx.obj = {"dask_cluster": client}  
        return ctx.obj

def main():
    cli_init.add_command(ingestion)
    cli_init.add_command(export)
    cli_init(obj={})

if __name__ == "__main__":
    main()
