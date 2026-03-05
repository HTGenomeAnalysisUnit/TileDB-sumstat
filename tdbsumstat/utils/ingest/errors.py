"""Custom exception types for the ingestion pipeline."""


class HarmonizationError(Exception):
    """Raised when data harmonization fails due to missing or invalid input."""
    pass
