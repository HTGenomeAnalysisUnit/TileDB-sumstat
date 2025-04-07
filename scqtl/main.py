#Create Dask cluster
from dask_jobqueue import LSFCluster as Cluster
from dask.distributed import LocalCluster
from dask.distributed import Client
import click
import cloup
from scqtl.cli.ingestion import ingestion
from scqtl.cli.export import export
@cloup.group(name="main", help="Single cell analysis with TileDB", no_args_is_help=True)
@cloup.option_group(
    "Dask cluster options",
    cloup.option("--workers", default = None, type=int, help = "Number of workers to use in the Dask cluster"),
    cloup.option("--cores_w", default = 4, type=int, help = "Number of core to use per single worker"),
    cloup.option("--memory_w", default = "8", type=str, help = "Number of GB to give as memory for each worker"),
	cloup.option("--local_cluster", default = False, is_flag = True, type=bool, help = "If you need to run on local cluster or distributed on different nodes"),

)

@click.pass_context
def cli_init(ctx, workers: int, cores_w: int, memory_w: str, local_cluster: bool):
	#Create Dask cluster
	if not workers is None:
		if local_cluster:
			cluster = LocalCluster(n_workers=workers, threads_per_worker=cores_w, memory_limit=f"{memory_w}GB", processes = 1)
		else:
			cluster = Cluster(cores = cores_w, memory=f"{memory_w}GB")
			cluster.scale(workers)
		client = Client(cluster)
		print(f"Cluster created: {cluster}")
		print(f"Dask dashboard available at {client.dashboard_link}")
		ctx.obj = {"dask_cluster": client}
		ctx.obj = {"workers": workers}
		return ctx.obj
	else:
		ctx.obj = {"workers": False}
		return ctx.obj

def main():
    cli_init.add_command(ingestion)
    cli_init.add_command(export)
    cli_init(obj={})

if __name__ == "__main__":
    main()
