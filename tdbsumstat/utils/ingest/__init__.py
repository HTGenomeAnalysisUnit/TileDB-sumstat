"""Ingestion package for TileDB-sumstat.

The ``Harmonize`` class is composed from focused mixin classes:

* :class:`~tdbsumstat.utils.ingest.schema.SchemaMixin`     – TileDB array creation
* :class:`~tdbsumstat.utils.ingest.mapping.MappingMixin`   – column-mapping CSV parsing
* :class:`~tdbsumstat.utils.ingest.harmonize.HarmonizeMixin` – data normalisation
* :class:`~tdbsumstat.utils.ingest.qc.QCMixin`             – optional gwaslab QC
* :class:`~tdbsumstat.utils.ingest.writer.WriterMixin`      – TileDB data writer
* :class:`~tdbsumstat.utils.ingest.metadata.MetadataMixin`  – metadata management
"""
import logging

from tdbsumstat.utils.ingest.errors import HarmonizationError
from tdbsumstat.utils.ingest.harmonize import HarmonizeMixin
from tdbsumstat.utils.ingest.mapping import MappingMixin
from tdbsumstat.utils.ingest.metadata import MetadataMixin
from tdbsumstat.utils.ingest.qc import QCMixin
from tdbsumstat.utils.ingest.schema import SchemaMixin
from tdbsumstat.utils.ingest.writer import WriterMixin

logger = logging.getLogger(__name__)

__all__ = ["Harmonize", "HarmonizationError"]


class Harmonize(SchemaMixin, MappingMixin, HarmonizeMixin, QCMixin, WriterMixin, MetadataMixin):
    """Pipeline class for harmonising and ingesting summary statistics into TileDB.

    Typical usage::

        h = Harmonize(mapping_file="mapping.csv", uri="my_tiledb",
                      type_sumstat="qtl", pvar_file=None,
                      type_trait="quant", mac=None, maf=None, permuted=False)

        # One-time setup
        h.create_tiledb()
        h.create_mapping()

        # Per-file loop
        for filepath, cell, gene in file_list:
            sumstat = pl.read_csv(filepath, separator="\\t", null_values="NA")
            h.harmonize(sumstat=sumstat, cell=cell, gene=gene, n=n, pheno_var=pheno_var)
            h.ingest_data(file_path=filepath)
            h.create_metadata(file_path=filepath)

        # Finalise
        h.merge_metadata_files()

    Parameters
    ----------
    mapping_file:
        Path to a header-less CSV mapping source column names to TileDB names.
    uri:
        Path where the TileDB array is (or will be) stored.
    type_sumstat:
        ``"gwas"`` or ``"qtl"``.
    pvar_file:
        Optional path to a pvar file used to align alleles.
    type_trait:
        ``"quant"`` or ``"binary"`` (used for GWAS only).
    mac:
        Minimum minor allele count filter applied during harmonisation.
    maf:
        Minimum minor allele frequency filter applied during harmonisation.
    permuted:
        If ``True``, SE is derived from the permuted p-value rather than
        being recalculated from BETA/SE.
    """

    def __init__(
        self,
        mapping_file: str,
        uri: str,
        type_sumstat: str,
        pvar_file: str,
        type_trait: str,
        mac: int,
        maf: float,
        permuted: bool,
    ) -> None:
        self.mapping_file = mapping_file
        self.uri = uri
        self.pvar_file = pvar_file
        self.type_sumstat = type_sumstat
        self.type_trait = type_trait
        self.mapping_types: dict = {}
        self.tiledb_types: dict = {}
        self.dimension_tiledb: list = []
        self.mac = mac
        self.maf = maf
        self.permuted = permuted
