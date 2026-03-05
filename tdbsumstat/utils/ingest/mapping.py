"""Column-mapping utilities for summary statistics ingestion."""
import logging

import pandas as pd

logger = logging.getLogger(__name__)


class MappingMixin:
    """Mixin that parses and validates the column-mapping CSV file."""

    def create_mapping(self) -> None:
        """Load the mapping file and populate ``self.mapping_types``.

        The mapping file is a header-less CSV with two columns:
        ``source_column_name, target_column_name``.

        Required target columns: ``BETA``, ``SE``.
        Required targets for position: either ``CHR`` + ``POS``, or ``SNPID``.

        Raises
        ------
        HarmonizationError
            If the file is empty, or required columns are missing.
        """
        from tdbsumstat.utils.ingest.errors import HarmonizationError

        df = pd.read_csv(self.mapping_file, header=None, names=["key", "value"])
        if df.empty:
            raise HarmonizationError("Mapping file is empty or not formatted correctly.")

        self.mapping_types = dict(zip(df["key"], df["value"]))

        if not all(col in self.mapping_types.values() for col in ["BETA", "SE"]):
            raise HarmonizationError("Mapping file must contain BETA and SE columns.")

        if not all(col in self.mapping_types.values() for col in ["CHR", "POS"]):
            if "SNPID" not in self.mapping_types.values():
                raise HarmonizationError("Mapping file must contain either CHR, and POS or SNPID columns.")

        logger.info("Column mapping loaded: %d column(s) mapped", len(self.mapping_types))
