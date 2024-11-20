
import tiledb
import numpy as np

def create_tiledb_schema(uri):
    pos_domain = (1, 3000000000)  # Example range for genomic positions
    dom = tiledb.Domain(
        tiledb.Dim(name="cell_type", dtype="ascii",  var=True, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]) ),
        tiledb.Dim(name="gene",dtype="ascii", var=True, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]) ),
        tiledb.Dim(name="position", domain = pos_domain, dtype=np.uint32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
    )
    schema = tiledb.ArraySchema(
        domain=dom,
        sparse=True,
        allows_duplicates=True,
        attrs=[
            tiledb.Attr(name="SNP", dtype="ascii", filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
	    tiledb.Attr(name="start_distance", dtype=np.int32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            tiledb.Attr(name="af", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            tiledb.Attr(name="beta", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            tiledb.Attr(name="se", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            tiledb.Attr(name="allele0", dtype="ascii", filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            tiledb.Attr(name="allele1", dtype="ascii", filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            tiledb.Attr(name="p-value", dtype=np.float64, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
        ]
    )
    tiledb.Array.create(uri, schema)


