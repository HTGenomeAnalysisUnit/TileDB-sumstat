"""Backward-compatibility shim.

All ingestion logic has been moved to the
:mod:`tdbsumstat.utils.ingest` package, which is composed of focused
modules:

* :mod:`tdbsumstat.utils.ingest.schema`    – TileDB schema creation
* :mod:`tdbsumstat.utils.ingest.mapping`   – column-mapping CSV parsing
* :mod:`tdbsumstat.utils.ingest.harmonize` – data normalisation
* :mod:`tdbsumstat.utils.ingest.qc`        – optional gwaslab QC
* :mod:`tdbsumstat.utils.ingest.writer`    – TileDB data writer
* :mod:`tdbsumstat.utils.ingest.metadata`  – metadata management

This file re-exports ``Harmonize`` and ``HarmonizationError`` so that
existing code that imports from ``tdbsumstat.utils.harmonize_ingest``
continues to work without modification.
"""
from tdbsumstat.utils.ingest import Harmonize, HarmonizationError  # noqa: F401

__all__ = ["Harmonize", "HarmonizationError"]
