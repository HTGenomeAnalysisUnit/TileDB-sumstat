import tiledb
import numpy as np
from .exceptions import HarmonizationError

def _create_tiledb(
                    type_sumstat:str,
                    uri:str
                    ):
        """Create the mapping and dtype definitions."""
        pos_domain = (1, 300000000)  # Example range for genomic positions
        chr_domain = (1, 24)  # Example range for genomic positions
        attrs=[
                tiledb.Attr(name="SNPID", dtype="ascii", filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="RSID", dtype="ascii", filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="EAF", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="BETA", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="SE", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Attr(name="P", dtype=np.float64, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
            ]
    
        if type_sumstat == "gwas":
            
            dimension_tiledb = ["CHR", "TRAIT", "POS"]
            dom = tiledb.Domain(
            tiledb.Dim(name="CHR", domain = chr_domain, dtype=np.uint16,  filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            tiledb.Dim(name="TRAIT", dtype="ascii",  var=True, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]) ),
            tiledb.Dim(name="POS", domain = pos_domain, dtype=np.uint32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
            )
        else:
            
            dimension_tiledb = ["CHR", "CELL", "GENE", "POS"]
            dom = tiledb.Domain(
                tiledb.Dim(name="CHR", domain = chr_domain, dtype=np.uint16,  filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
                tiledb.Dim(name="CELL", dtype="ascii",  var=True, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]) ),
                tiledb.Dim(name="GENE",dtype="ascii", var=True, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]) ),
                tiledb.Dim(name="POS", domain = pos_domain, dtype=np.uint32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
                )
            
            attrs = attrs + [
            tiledb.Attr(name="DIST", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]))
            ]
        schema = tiledb.ArraySchema(
                domain=dom,
                attrs=attrs,
                sparse=True,
                allows_duplicates=False
        )
        tiledb.Array.create(uri, schema)