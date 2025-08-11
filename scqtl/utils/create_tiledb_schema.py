import tiledb
import numpy as np

def create_tiledb_schema(uri:str, type_sumstat:str):
    pos_domain = (1, 300000000)  # Example range for genomic positions
    chr_domain = (1, 22)  # Example range for genomic positions
    #define the domain for chromosome, cell, gene and position
    if type_sumstat not in ["gwas", "scqtl"]:
        raise ValueError("type_sumstat must be either 'gwas' or 'scqtl'")
    if type_sumstat == "gwas":
        dom = tiledb.Domain(
            tiledb.Dim(name="CHR", domain=chr_domain, dtype=np.uint16, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            tiledb.Dim(name="TRAIT", dtype="ascii", var=True, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            tiledb.Dim(name="POS", domain=pos_domain, dtype=np.uint32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
        )
    elif type_sumstat == "scqtl":
        # For scQTL, we need to define additional dimensions for CELL and GENE
        dom = tiledb.Domain(
                tiledb.Dim(name="CHR", domain=chr_domain, dtype=np.uint16, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Dim(name="CELL", dtype="ascii", var=True, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Dim(name="GENE", dtype="ascii", var=True, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Dim(name="POS", domain=pos_domain, dtype=np.uint32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
            )
    attributes=[
            tiledb.Attr(name="SNP", dtype="ascii", filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),	        
            tiledb.Attr(name="AF", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            tiledb.Attr(name="BETA", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            tiledb.Attr(name="SE", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            tiledb.Attr(name="P", dtype=np.float64, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            tiledb.Attr(name="N", dtype=np.int64, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
        ]
    if type_sumstat == "scqtl":
         attributes = attributes + [
            tiledb.Attr(name="DIST", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
            ]
    schema = tiledb.ArraySchema(
        domain=dom,
        sparse=True,
        allows_duplicates=True,
        attrs = attributes    
    )
    tiledb.Array.create(uri, schema)


