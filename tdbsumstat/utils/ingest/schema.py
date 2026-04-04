"""TileDB array schema creation for summary statistics ingestion."""
import logging

import numpy as np
import tiledb

logger = logging.getLogger(__name__)


class SchemaMixin:
    """Mixin that provides TileDB array schema creation."""

    def create_tiledb(self) -> None:
        """Create the TileDB array with the appropriate schema.

        The schema varies by ``self.type_sumstat``:

        * ``"gwas"`` – dimensions: CHR, TRAIT, POS
        * ``"qtl"``  – dimensions: CHR, CELL, GENE, POS  (+ DIST attribute)
        """
        pos_domain = (1, 300_000_000)
        chr_domain = (1, 24)

        attrs = [
            tiledb.Attr(name="SNPID", dtype="ascii", filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            tiledb.Attr(name="RSID", dtype="ascii", filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            tiledb.Attr(name="EAF", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            tiledb.Attr(name="BETA", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            tiledb.Attr(name="SE", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            tiledb.Attr(name="P", dtype=np.float64, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
        ]

        if self.type_sumstat == "gwas":
            self.dimension_tiledb = ["CHR", "TRAIT", "POS"]
            dom = tiledb.Domain(
                tiledb.Dim(
                    name="CHR", domain=chr_domain, dtype=np.uint16,
                    filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]),
                ),
                tiledb.Dim(
                    name="TRAIT", dtype="ascii", var=True,
                    filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]),
                ),
                tiledb.Dim(
                    name="POS", domain=pos_domain, dtype=np.uint32,
                    filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]),
                ),
            )
        else:
            self.dimension_tiledb = ["CHR", "CELL", "GENE", "POS"]
            dom = tiledb.Domain(
                tiledb.Dim(
                    name="CHR", domain=chr_domain, dtype=np.uint16,
                    filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]),
                ),
                tiledb.Dim(
                    name="CELL", dtype="ascii", var=True,
                    filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]),
                ),
                tiledb.Dim(
                    name="GENE", dtype="ascii", var=True,
                    filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]),
                ),
                tiledb.Dim(
                    name="POS", domain=pos_domain, dtype=np.uint32,
                    filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)]),
                ),
            )
            attrs = attrs + [
                tiledb.Attr(name="DIST", dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=5)])),
            ]

        schema = tiledb.ArraySchema(
            domain=dom,
            attrs=attrs,
            sparse=True,
            allows_duplicates=False,
        )
        tiledb.Array.create(self.uri, schema)
        logger.info("TileDB array created at %s (type_sumstat=%s)", self.uri, self.type_sumstat)
