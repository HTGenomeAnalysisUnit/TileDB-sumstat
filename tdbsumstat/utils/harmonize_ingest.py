from tdbsumstat.harmonize import _create_mapping, _align_alleles, _fill_missing_columns,_qc_sumstat
from tdbsumstat.ingest import _create_metadata, _create_tiledb, _export_metadata_to_csv,_ingest_data,_merge_metadata_files, _write_individual_metadata
from pathlib import Path
class Harmonize:
    def __init__(
            self,
            mapping_file: Path,
            uri: str,
            type_sumstat: str,
            type_trait: str,
            mac: int,
            pvar_file: Path | None = None):
        mapping_file = mapping_file
        uri = uri
        pvar_file = pvar_file
        type_sumstat = type_sumstat
        type_trait = type_trait
        mac = mac
      
