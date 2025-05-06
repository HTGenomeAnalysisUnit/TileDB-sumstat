import pandas as pd

existing_gene_celltype = {}
        with tiledb.open(uri, mode="r") as A:
            if "GENE_CELLTYPE" in A.meta:
                existing_gene_celltype = json.loads(A.meta["GENE_CELLTYPE"])

        with tiledb.open(uri, mode="w") as A:
            if not existing_gene_celltype:  # If GENE_CELLTYPE doesn't exist yet
                new_dict = {file[1]: chunk_pl["phenotype_id"].unique().to_list()}
                A.meta["GENE_CELLTYPE"] = json.dumps(new_dict)
            else:
                if file[1] not in existing_gene_celltype:
                    existing_gene_celltype[file[1]] = chunk_pl["phenotype_id"].unique().to_list()
                else:
                    # Append new cell types to the existing list
                    existing_gene_celltype[file[1]].extend(chunk_pl["phenotype_id"].unique().to_list())
                A.meta["GENE_CELLTYPE"] = json.dumps(existing_gene_celltype)
