from .exceptions import HarmonizationError
import pandas as pd
from pathlib import Path

def _create_mapping(
          sumstat: pd.DataFrame,
          mapping_file: Path | None = None
          ):
        if mapping_file is None:
             raise HarmonizationError("No mapping file provided. Provide a mapping file")
        if not mapping_file.is_file():
            raise FileNotFoundError(f"Mapping file {mapping_file} does not exist")
        
        df = pd.read_csv(mapping_file, header=None, names=["key", "value"])
        if df.empty:
            raise HarmonizationError("Mapping file is empty or not formatted correctly.")
        mapping_types = dict(zip(df["key"], df["value"]))
        #check that "BETA", "SE" are in the vlaues of the mapping_types
        if not all(col in mapping_types.values() for col in ["BETA", "SE"]):
            raise HarmonizationError("Mapping file must contain BETA and SE columns.")
        # Check if CHR and POS or SNP are present
        if not all(col for col in ["CHR", "POS"] if col in mapping_types.values()):
            if "SNPID" not in mapping_types.values():
                raise HarmonizationError("Mapping file must contain either CHR, and POS or SNPID columns.")
        chunk_pl = sumstat.rename(mapping_types)
        return chunk_pl