import os
import shutil
import tempfile

import pybedtools


class bed_parser:
    def __init__(
        self,
        mod_code=None,
        min_valid_cov=0,
        bedgraph=False,
        sat_type=None,
        pre_subset_censat=False,
        temp_dir=None,
        cache=True,
    ):
        """
        Initialize the parser with optional filtering parameters

        Args:
            min_valid_cov (int): Minimum coverage threshold
            mod_code (str): Modification code to filter
            sat_type (list): Satellite types to filter
            temp_dir (str): Directory for temporary file storage
            cache (bool): Enable caching of BEDTools operations
        """
        self.mod_code = mod_code
        self.min_valid_cov = min_valid_cov
        self.bedgraph = bedgraph
        self.sat_type = sat_type or []
        self.pre_subset_censat = pre_subset_censat

        # Set up temporary directory for intermediate files
        self.temp_dir = temp_dir or tempfile.gettempdir()

        # Enable caching to improve performance for repeated operations
        if cache:
            pybedtools.set_tempdir(self.temp_dir)

    def check_bedtools_installed(self):
        if shutil.which("bedtools") is None:
            raise EnvironmentError(
                "bedtools is required but not installed. Please install it before using this package."
            )

    def read_and_filter_bedfile(self, filepath, file_type="bedmethyl"):
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
        if file_type == "bedmethyl":
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

                def bedmethyl_filter(feature):
                    try:
                        if (
                            feature[3] == self.mod_code
                            and float(feature[4]) >= self.min_valid_cov
                        ):
                            # Assuming 4 columns, last column is methylation value
                            return feature
                    except (IndexError, ValueError):
                        return False

                    return False  # Explicitly return False for features not meeting criteria

                bedtool = bedtool.filter(bedmethyl_filter).cut([0, 1, 2, 10])

        elif file_type == "censat":
            if self.pre_subset_censat:
                return bedtool
            else:
                # Filter by satellite type
                def censat_filter(feature):
                    try:
                        return any(
                            sat.lower() in feature[3].lower() for sat in self.sat_type
                        )
                    except (IndexError, AttributeError):
                        return False

                bedtool = bedtool.filter(censat_filter)

        return bedtool.saveas()

    def bedtool_to_chrom_dict(self, bedtool):
        """
        Group intersected intervals by chromosome

        Args:
            intersected_bedtool (pybedtools.BedTool): Intersected intervals

        Returns:
            dict: Chromosome-wise dictionary of methylation data
        """
        # Convert to DataFrame for easier manipulation
        df = bedtool.to_dataframe()

        # Group by chromosome
        return {chrom: group for chrom, group in df.groupby(df.columns[0])}

    def process_files(self, bedmethyl_path, censat_path):
        """
        Main processing method to intersect and group files

        Args:
            bedmethyl_path (str): Path to BedMethyl file
            censat_path (str): Path to CenSat file

        Returns:
            dict: Chromosome-wise dictionary of methylation data
        """
        # Read and filter files
        bedmethyl_bedtool = self.read_and_filter_bedfile(bedmethyl_path, "bedmethyl")
        censat_bedtool = self.read_and_filter_bedfile(censat_path, "censat")

        # Perform intersection
        # -wa: Write original A entry
        intersected_bedmethyl_bedtool = bedmethyl_bedtool.intersect(
            censat_bedtool, wa=True
        )

        return self.bedtool_to_chrom_dict(
            intersected_bedmethyl_bedtool
        ), self.bedtool_to_chrom_dict(censat_bedtool)
