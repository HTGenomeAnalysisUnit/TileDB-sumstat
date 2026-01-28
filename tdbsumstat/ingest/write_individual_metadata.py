import json
import os
from pathlib import Path
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def _write_individual_metadata(
            metadata,
            file_path,
            uri,):
        """Write metadata to individual JSON file."""
        metadata_dir = f"{uri}_metadata_parts"
        os.makedirs(metadata_dir, exist_ok=True)
    
        # Create a safe filename from the original file path
        file_stem = Path(file_path).stem
        metadata_file = os.path.join(metadata_dir, f"{file_stem}.json")
    
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
    
        logger.info(f"Metadata written to {metadata_file}")