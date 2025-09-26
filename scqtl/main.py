#Create Dask cluster
from dask_jobqueue import LSFCluster
from dask.distributed import LocalCluster
from dask.distributed import Client
import click
import cloup
import os
from scqtl.cli.ingestion import ingest
from scqtl.cli.export import export

@cloup.group(name="scqtl", help="Single cell analysis with TileDB", no_args_is_help=True)
@cloup.option_group(
    "Dask cluster options",
    cloup.option("--workers", default = None, type=int, help = "Number of workers to use in the Dask cluster"),
    cloup.option("--cores_w", default = 4, type=int, help = "Number of core to use per single worker"),
    cloup.option("--memory_w", default = "12", type=str, help = "Number of GB to give as memory for each worker"),
	cloup.option("--type_cluster", default = "local",  help = "If you need to run on local cluster or distributed on different nodes"),

)

@click.pass_context
def cli_init(ctx, workers: int, cores_w: int, memory_w: str, type_cluster: str):
	#Create Dask cluster
	if not workers is None:
		if type_cluster == 'local':
			cluster = LocalCluster(n_workers=workers, threads_per_worker=cores_w, memory_limit=f"{memory_w}GB", processes = True, dashboard_address=':0')
		else:
			cluster = LSFCluster(
        		queue='normal',
        		walltime='00:40',
        		log_directory='{}/dask_logs'.format(os.getcwd()),
        		cores=2,
				processes = 1,
        		memory='{} Gb'.format(20),
        		mem=20*1e+9,  # should be in bytes
        		lsf_units='mb',
        		job_extra=[
            		'-G team151',
            		'-g /lt9/dask',
            		'-R "select[mem>{}] rusage[mem={}]"'.format(
                		int(20*1e+3),
                		int(20*1e+3)
            		)
        			],
        		use_stdin=True
    			)
			#cluster = LSFCluster(cores = cores_w, memory=f"{memory_w}GB")
			cluster.scale(workers)
		client = Client(cluster)
		print(f"Cluster created: {cluster}")
		print(f"Dask dashboard available at {client.dashboard_link}")
		ctx.obj = {
    		"dask_cluster": client, 
    		"workers": workers,
    		"cluster_object": cluster  # Also store the cluster itself
			}
		return ctx.obj
	else:
		ctx.obj = {"workers": False, "dask_cluster": None}
		return ctx.obj

def main():
    cli_init.add_command(ingest)
    cli_init.add_command(export)
    cli_init(obj={})

if __name__ == "__main__":
    main()
