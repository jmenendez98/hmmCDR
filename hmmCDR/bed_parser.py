import pandas as pd
import numpy as np
import pybedtools
import os
import tempfile
import shutil

class bed_parser:
    def __init__(self, min_valid_cov=0, mod_code=None, bedgraph=False, sat_type=None, 
                 temp_dir=None, cache=True):
        """
        Initialize the parser with optional filtering parameters
        
        Args:
            min_valid_cov (int): Minimum coverage threshold
            mod_code (str): Modification code to filter
            sat_type (list): Satellite types to filter
            temp_dir (str): Directory for temporary file storage
            cache (bool): Enable caching of BEDTools operations
        """
        self.min_valid_cov = min_valid_cov
        self.mod_code = mod_code
        self.sat_type = sat_type or []
        self.bedgraph = bedgraph
        
        # Set up temporary directory for intermediate files
        self.temp_dir = temp_dir or tempfile.gettempdir()
        
        # Enable caching to improve performance for repeated operations
        if cache:
            pybedtools.set_tempdir(self.temp_dir)

    def check_bedtools_installed(self):
        if shutil.which("bedtools") is None:
            raise EnvironmentError("bedtools is required but not installed. Please install it before using this package.")

    def read_and_filter_bedfile(self, filepath, file_type='bedmethyl'):
        """
        Read and filter BED-like files with advanced optimizations
        
        Args:
            filepath (str): Path to the input file
            file_type (str): Type of file ('bedmethyl' or 'censat')
        
        Returns:
            pybedtools.BedTool: Filtered genomic intervals
        """
        # Validate file existence
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")
        
        # Create BedTool directly from file
        bedtool = pybedtools.BedTool(filepath)

        # Filtering logic based on file type
        if file_type == 'bedmethyl':
            if self.bedgraph:
                # Filter by bedgraph code and coverage
                def bedgraph_filter(feature):
                    try:
                        # Assuming 4 columns, last column is methylation value
                        return feature
                    except (IndexError, ValueError):
                        return False
                bedtool = bedtool.filter(bedgraph_filter)

            else:
                # Original BedMethyl filtering logic
                def bedmethyl_filter(feature):
                    try:
                        # Check filtering conditions
                        return (
                            feature[3] == self.mod_code
                            and float(feature[4]) >= self.min_valid_cov
                        )
                    except (IndexError, ValueError):
                        return False
                bedtool = bedtool.filter(bedmethyl_filter)

        elif file_type == 'censat' and self.sat_type:
            # Filter by satellite type
            def censat_filter(feature):
                try:
                    return any(sat.lower() in feature[3].lower() 
                               for sat in self.sat_type)
                except (IndexError, AttributeError):
                    return False
            bedtool = bedtool.filter(censat_filter)
        
        return bedtool
    
    def intersect_files(self, bedmethyl_path, censat_path):
        """
        Perform intersection between BedMethyl and CenSat files
        
        Args:
            bedmethyl_path (str): Path to BedMethyl file
            censat_path (str): Path to CenSat file
        
        Returns:
            pybedtools.BedTool: Intersected genomic intervals
            pybedtools.BedTool: CenSat intervals used for intersection
        """
        # Read and filter files
        bedmethyl_tool = self.read_and_filter_bedfile(bedmethyl_path, 'bedmethyl')
        censat_tool = self.read_and_filter_bedfile(censat_path, 'censat')

        # Perform intersection
        # -u: Report original A entry once if overlaps with any B
        # -wa: Write original A entry
        intersected = bedmethyl_tool.intersect(censat_tool, u=True, wa=True)

        return intersected

    def group_by_chromosome(self, methylation_bedtool):
        """
        Group intersected intervals by chromosome
        
        Args:
            intersected_bedtool (pybedtools.BedTool): Intersected intervals
        
        Returns:
            dict: Chromosome-wise dictionary of methylation data
        """
        # Convert to DataFrame for easier manipulation
        if self.bedgraph:
            methyl_df = methylation_bedtool.to_dataframe(names=np.arange(0,4,1))
        else:
            methyl_df = methylation_bedtool.to_dataframe(names=np.arange(0,18,1))
            methyl_df = methyl_df.iloc[:, [0, 1, 2, 10]]

        if methyl_df.empty:
            raise ValueError("The intersection resulted in an empty DataFrame.")

        # Group by chromosome
        return {chrom: group for chrom, group in methyl_df.groupby(methyl_df.columns[0])}
    
    def process_files(self, bedmethyl_path, censat_path):
        """
        Main processing method to intersect and group files
        
        Args:
            bedmethyl_path (str): Path to BedMethyl file
            censat_path (str): Path to CenSat file
        
        Returns:
            dict: Chromosome-wise dictionary of methylation data
        """
        # Intersect files
        methylation_bedtool = self.intersect_files(bedmethyl_path, censat_path)

        return self.group_by_chromosome(methylation_bedtool)