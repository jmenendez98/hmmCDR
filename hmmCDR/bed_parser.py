import os
import numpy as np
from typing import Dict, Optional, List, Union
import tempfile

class bed_parser:
    """hmmCDR parser to read in region and methylation bed files."""

    def __init__(
        self,
        mod_code: Optional[str] = None,
        min_valid_cov: int = 0,
        methyl_bedgraph: bool = False,
        sat_type: Optional[Union[str, List[str]]] = None,
        edge_filter: int = 50000,
        regions_prefiltered: bool = False
    ):
        """
        Initialize the parser with optional filtering parameters.

        Args:
            mod_code: Modification code to filter
            min_valid_cov: Minimum coverage threshold
            methyl_bedgraph: Whether the file is a bedgraph
            sat_type: Satellite type(s) to filter
            edge_filter: Amount to remove from edges of active_hor regions
            regions_prefiltered: Whether the regions bed is already subset
        """
        self.mod_code = mod_code
        self.min_valid_cov = min_valid_cov
        self.methyl_bedgraph = methyl_bedgraph
        self.sat_type = [sat_type] if isinstance(sat_type, str) else (sat_type or [])
        self.edge_filter = edge_filter
        self.regions_prefiltered = regions_prefiltered
        self.temp_dir = tempfile.gettempdir()

    def read_and_filter_regions(self, regions_path: str) -> Dict[str, Dict[str, np.ndarray]]:
        """
        Read and filter regions from a BED file.
        
        Args:
            regions_path: Path to the regions BED file
            
        Returns:
            Dictionary mapping chromosomes to their start/end positions
            
        Raises:
            FileNotFoundError: If regions_path doesn't exist
            TypeError: If BED file is incorrectly formatted
        """
        if not os.path.exists(regions_path):
            raise FileNotFoundError(f"File not found: {regions_path}")

        region_dict: Dict[str, Dict[str, np.ndarray]] = {}

        with open(regions_path, 'r') as file:
            lines = [line.strip().split('\t') for line in file]
            
            if any(len(cols) < 3 for cols in lines):
                raise TypeError(f"Less than 3 columns in {regions_path}. Likely incorrectly formatted bed file.")

            chrom_lines = {}
            for cols in lines:
                if ((self.regions_prefiltered) or (len(cols) > 3 and cols[3] in self.sat_type)):
                    chrom = cols[0]
                    if chrom not in chrom_lines:
                        chrom_lines[chrom] = []
                    chrom_lines[chrom].append(cols)

            for chrom, chrom_data in chrom_lines.items():
                starts = np.array([int(cols[1]) + self.edge_filter for cols in chrom_data], dtype=int)
                ends = np.array([int(cols[2]) - self.edge_filter for cols in chrom_data], dtype=int)
                region_dict[chrom] = {"starts": starts, "ends": ends}
                    
        return region_dict
    
    def read_and_filter_methylation(self, methylation_path):
        """
        Read and filter methylation data from a BED file.
        
        Args:
            methylation_path: Path to the methylation BED file
            
        Returns:
            Dictionary mapping chromosomes to their methylation data
            
        Raises:
            FileNotFoundError: If methylation_path doesn't exist
            TypeError: If BED file is incorrectly formatted
            ValueError: If trying to filter bedgraph by coverage
        """
        if not os.path.exists(methylation_path):
            raise FileNotFoundError(f"File not found: {methylation_path}")
        
        if self.methyl_bedgraph and self.min_valid_cov > 0:
            raise ValueError(f"{methylation_path} bedgraph file cannot be filtered by coverage.")

        methylation_dict: Dict[str, Dict[str, np.ndarray]] = {}

        with open(methylation_path, 'r') as file:
            lines = [line.strip().split('\t') for line in file]
            
            if any(len(cols) < (4 if self.methyl_bedgraph else 11) for cols in lines):
                raise TypeError(f"Insufficient columns in {methylation_path}. Likely incorrectly formatted.")

            chrom_lines = {}
            for cols in lines:
                chrom = cols[0]
                if self.methyl_bedgraph or (cols[3] == self.mod_code and float(cols[4]) >= self.min_valid_cov):
                    if chrom not in chrom_lines:
                        chrom_lines[chrom] = []
                    chrom_lines[chrom].append(cols)

            for chrom, chrom_data in chrom_lines.items():
                starts = np.array([int(cols[1]) for cols in chrom_data], dtype=int)
                scores = np.array([float(cols[3 if self.methyl_bedgraph else 10]) for cols in chrom_data], dtype=float)
                methylation_dict[chrom] = {"starts": starts, "scores": scores}
                    
        return methylation_dict

    def process_files(self, methylation_path, regions_path):
        """
        Process and intersect methylation and regions files.
        
        Args:
            methylation_path: Path to methylation BED file
            regions_path: Path to regions BED file
            
        Returns:
            Tuple of (region_dict, filtered_methylation_dict)
        """
        region_dict = self.read_and_filter_regions(regions_path)
        methylation_dict = self.read_and_filter_methylation(methylation_path)

        filtered_methylation_dict = {}

        for chrom, regions in region_dict.items():
            if chrom not in methylation_dict:
                continue

            methylation_data = methylation_dict[chrom]
            
            # Vectorized overlap check
            region_starts = regions['starts'][:, np.newaxis]
            region_ends = regions['ends'][:, np.newaxis]
            methyl_starts = methylation_data['starts']
            methyl_ends = methyl_starts + 1

            # Check overlap
            overlaps = (methyl_starts <= region_ends) & (methyl_ends >= region_starts)
            valid_positions = np.any(overlaps, axis=0)

            if np.any(valid_positions):
                filtered_methylation_dict[chrom] = {
                    "starts": methylation_data['starts'][valid_positions],
                    "scores": methylation_data['scores'][valid_positions]
                }

        return region_dict, filtered_methylation_dict