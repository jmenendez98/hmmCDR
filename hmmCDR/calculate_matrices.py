import pandas as pd
import numpy as np 
import pybedtools
import argparse
import os
import concurrent.futures

from hmmCDR.bed_parser import bed_parser

class calculate_matrices:
    def __init__(self, window_size, step_size, min_prior_size,
                 enrichment, output_label):

        self.window_size = window_size
        self.step_size = step_size
        self.min_prior_size = min_prior_size

        self.enrichment = enrichment
        self.output_label = output_label
        
        self.MAX_RETRIES = 10
        

    def create_chrom_windows(self, methylation_bedtool, gap_threshold=10000):
        """
        Create fixed-size windows across chromosomes, separating windows in regions with large gaps.

        Args:
            methylation_bedtool (pybedtools.BedTool): Input bedtool with 'chrom', 'start', and 'end' columns.
            gap_threshold (int): Maximum allowable gap between intervals before splitting into separate regions.

        Returns:
            pybedtools.BedTool: BedTool object containing fixed-size windows.
        """

        # Merge intervals that are closer than the gap threshold
        merged_bedtool = methylation_bedtool.merge(d=gap_threshold)

        # Generate fixed-size windows from the merged intervals
        windows_bedtool = merged_bedtool.window_maker(b=merged_bedtool, w=self.window_size, s=self.step_size)

        return windows_bedtool

    def mean_within_windows(self, methylation_bedtool, windows_bedtool):
        """
        Calculates mean methylation values within each of the windows

        Args:
            methylation_bedgraph (pd.DataFrame): Input DataFrame with 'chrom', 'start', and 'end' columns.
            windows_bedtool (pybedtools.BedTool): BedTool object containing fixed-size windows.

        Returns:
            pybedtools.BedTool: BedTool object containing fixed-size windows with methylation as the fourth column.
        """
        # Find 
        windows_means_bedtool = windows_bedtool.map(methylation_bedtool, c=4, o='mean')
        return windows_means_bedtool

    def find_prior_percentile(self, windows_means_bedtool, prior_threshold):
        """
        Calculates threshold if percentiles are being used in generating priors

        Args:
            windows_means_bedtool (pybedtools.BedTool): BedTool object containing fixed-size windows and their mean methylation.
            windows_means_bedtool (float): Input float to use as percentile to find. 

        Returns:
            float: Threshold to be used to determine if prior is a valid CDR.
        """
        percentile = prior_threshold
        window_mean_values = pd.to_numeric(windows_means_bedtool.to_dataframe().iloc[:, 3].replace('.', np.nan)).dropna()
        return np.percentile(window_mean_values, q=percentile)

    def find_priors(self, windows_means_bedtool, filter_threshold):
        """
        Create a DataFrame of genomic regions filtered by mean values and merged based on conditions.

        Args:
            windows_means_bedtool (pd.pybedtools.BedTool): BedTool object containing fixed-size windows and their mean methylation.
            filter_threshold (float): Threshold for filtering windows.
            
        Returns:
            pd.DataFrame: DataFrame of filtered and merged genomic windows.
        """
        # Filter windows based on enrichment or depletion
        def safe_filter(feature):
            try:
                # Skip entries with '.' or empty values
                if feature[3] == '.' or feature[3] == '':
                    return False
                
                # Convert to float and apply filtering
                value = float(feature[3])
                
                if self.enrichment:
                    return value > filter_threshold
                else:
                    return value < filter_threshold
            
            except (ValueError, IndexError, TypeError):
                # Handle any conversion or index errors
                return False
            
        filtered_windows_bedtool = windows_means_bedtool.filter(safe_filter)

        # Merge adjacent or overlapping regions
        merged_windows_bedtool = filtered_windows_bedtool.merge()

        # Filter by minimum size using pybedtools filter
        final_windows_bedtool = merged_windows_bedtool.filter(
            lambda x: int(x.end) - int(x.start) >= self.min_prior_size
        )

        return final_windows_bedtool
    
    def calculate_emission_thresholds(self, bed4Methyl):
        if not self.hmm_percentile_emissions:
            return sorted([self.w, self.x, self.y, self.z])
        
        methylation_scores = pd.to_numeric(bed4Methyl['name'].replace('.', np.nan), errors='coerce').dropna()
        methylation_scores = [0] + methylation_scores[methylation_scores != 0].tolist()
        
        return sorted(np.percentile(methylation_scores, q=[self.w, self.x, self.y, self.z]))

    def assign_emissions(self, bed4Methyl, emission_thresholds):
        """
        Assign emission states using vectorized operations.
        
        Args:
            bed4Methyl (pd.DataFrame): DataFrame with methylation values
            emission_thresholds (list): Sorted list of threshold values
            
        Returns:
            pd.DataFrame: DataFrame with added 'emission' column
        """
        # Convert values to numeric array
        values = pd.to_numeric(bed4Methyl['name'], errors='coerce')
        
        # Create a matrix of boolean comparisons
        comparison_matrix = values.reshape(-1, 1) <= np.array(emission_thresholds)
        
        # Get the first True index for each row (vectorized operation)
        bed4Methyl['emission'] = np.argmax(comparison_matrix, axis=1)
        
        # Handle cases where all comparisons are False
        max_emission = len(emission_thresholds) - 1
        bed4Methyl.loc[~comparison_matrix.any(axis=1), 'emission'] = max_emission
        
        return bed4Methyl
    
    def assign_priors(self, methylation_bedtool, cdr_prior_bedtool):
        """
        Create a DataFrame of genomic regions filtered by mean values and merged based on conditions.

        Args:
            methylation_bedtool (pd.pybedtools.BedTool): BedTool object containing CpG positions and their methylation.
            cdr_priors (pd.pybedtools.BedTool): BedTool object containing regions to be used as priors.
            
        Returns:
            pd.Dataframe: Pandas dataframe chrom positions, fraction modified, and prior state.
        """
        methylation_in_priors_bedtool = methylation_bedtool.intersect(
            cdr_prior_bedtool, 
            wao=True,  # Write all original entries
        )

        # Convert methylation bedtool with priors assigned to pd.Dataframe
        methylation_w_priors_df = methylation_in_priors_bedtool.to_dataframe(names=np.arange(0,8,1))
        methylation_w_priors_df = methylation_w_priors_df.iloc[:, [0,1,2,3,7]]
        methylation_w_priors_df.columns = ['chrom', 'start', 'end' ,'fractionmodified', 'prior']

        return methylation_w_priors_df

    def calculate_transition_matrix(self, methylation_w_priors_df):
        """
        Calculate transition matrix for binary states using vectorized operations.
        
        Args:
            methylation_w_priors_df (pd.DataFrame): DataFrame with 'prior' column containing binary states
            
        Returns:
            np.ndarray: 2x2 transition probability matrix
        """
        # Convert priors to numpy array for faster operations
        states = methylation_w_priors_df.iloc[: , 4].to_numpy()
        
        # Create pairs of consecutive states using array slicing
        state_pairs = np.vstack((states[:-1], states[1:])).T
    
        # Use numpy's unique with return_counts to get transition counts
        unique_pairs, counts = np.unique(state_pairs, axis=0, return_counts=True)
    
        # Initialize transition matrix with zeros
        transition_matrix = np.zeros((2, 2))
        
        # Fill transition matrix using unique pairs and counts
        for (i, j), count in zip(unique_pairs, counts):
            transition_matrix[int(i), int(j)] = count
        
        # Normalize by row sums (avoiding division by zero)
        row_sums = transition_matrix.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1  # Avoid division by zero
        transition_matrix = transition_matrix / row_sums

        return transition_matrix

    def priors_single_chromosome(self, chrom, chrom_methylation_df, percentile, threshold):
        """
        Find prior CDRs for a single chromosome with adaptive thresholding.

        Args:
            chrom (str): Chromosome name
            chrom_methylation_df (pd.DataFrame): Methylation data for the chromosome
            percentile (bool): Whether to use percentile-based thresholding
            threshold (float): Initial threshold or percentile value

        Returns:
            tuple: Chromosome name, prior CDRs DataFrame, window mean BedTool
        """
        methylation_bedtool = pybedtools.BedTool.from_dataframe(chrom_methylation_df) # Convert the input dataFrame to a bedTool object
        window_bedtool = self.create_chrom_windows(methylation_bedtool) # Create windows
        window_mean_bedtool = self.mean_within_windows(methylation_bedtool, window_bedtool) # Find mean within windows

        for attempt in range(self.MAX_RETRIES):
            # Determine current CDR threshold
            if percentile: 
                prior_threshold = self.find_prior_percentile(window_mean_bedtool, threshold)
            else:
                prior_threshold = threshold

            # Find priors using current threshold
            cdr_prior_bedtool = self.find_priors(window_mean_bedtool, prior_threshold)

            # assign emissions to methylation bedtool
            methylation_w_emission_df = self.assign_emissions(methylation_bedtool, emission_thresholds)

            # add priors on to methylation bedtool
            methylation_w_priors_df = self.assign_priors(methylation_bedtool, cdr_prior_bedtool)

            # calculate emission and transition matrices with assigned priors
            emission_matrix = self.calculate_emission_matrix(methylation_w_priors_df)
            transition_matrix = self.calculate_transition_matrix(methylation_w_priors_df)

            return chrom, cdr_prior_bedtool, window_mean_bedtool
            
            # If no priors found, adjust threshold
            print(f'No prior CDR Found on {chrom} with threshold - {threshold}')
        
        # If we've exhausted all retry attempts
        raise RuntimeError(
            f"Failed to detect prior subCDRs for {chrom} after {self.MAX_RETRIES} retries."
            f"with final percentile - {prior_threshold}"
        )

    def priors_all_chromosomes(self, bed4Methyl_chrom_dict, percentile, threshold):
        hmmCDRpriors_chrom_dict = {}
        windowsmean_chrom_dict = {}
        chromosomes = bed4Methyl_chrom_dict.keys()

        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = {
                executor.submit(
                    self.priors_single_chromosome, 
                    chrom,
                    bed4Methyl_chrom_dict[chrom], 
                    percentile, 
                    threshold
                ): chrom for chrom in chromosomes
            }

            for future in concurrent.futures.as_completed(futures):
                chrom, priors_df, windows_mean_df = future.result()
                hmmCDRpriors_chrom_dict[chrom] = priors_df
                windowsmean_chrom_dict[chrom] = windows_mean_df

        self.chromosomes = chromosomes
        self.hmmCDRpriors_chrom_dict = hmmCDRpriors_chrom_dict
        self.windowsmean_chrom_dict = hmmCDRpriors_chrom_dict

        return hmmCDRpriors_chrom_dict, windowsmean_chrom_dict


def main():
    argparser = argparse.ArgumentParser(description='Process bedMethyl and CenSat BED file to produce hmmCDR priors')
    
    # Required arguments
    argparser.add_argument('bedmethyl_path', type=str, help='Path to the bedMethyl file')
    argparser.add_argument('censat_path', type=str, help='Path to the CenSat BED file')
    argparser.add_argument('output_path', type=str, help='Path to the output priorCDRs BED file')
    
    # Parser Arguments
    argparser.add_argument('-m', '--mod_code', type=str, default='m', help='Modification code to filter bedMethyl file (default: "m")')
    argparser.add_argument('-s', '--sat_type', type=str, default='H1L', help='Comma-separated list of satellite types/names to filter CenSat bed file. (default: "H1L")')
    argparser.add_argument('--bedgraph', action='store_true', default=False, help='Flag indicating if the input is a bedgraph. (default: False)')
    argparser.add_argument('--min_valid_cov', type=int, default=10, help='Minimum valid coverage to consider a methylation site(read from full modkit pileup files). (default: 10)')
    
    # calculate_matrices Arguments
    argparser.add_argument('--window_size', type=int, default=1190, help='Window size to calculate prior regions. (default: 1190)')
    argparser.add_argument('--step_size', type=int, default=1190, help='Step size when calculation windows for priors. (default: 1190)')
    argparser.add_argument('--prior_threshold', type=float, default=30.0, help='Threshold for determining if a window is a CDR. Uses this percentile if -prior_use_percentile is passed (default: 30.0)')
    argparser.add_argument('--prior_use_percentile', action='store_true', default=False, help='Whether or not to use percentile when calculating priors. (default: False)')
    argparser.add_argument('--min_prior_size', type=int, default=8330, help='Minimum size for CDR regions. (default: 8330)')
    argparser.add_argument('--enrichment', action='store_true', default=False, help='Enrichment flag. Pass in if you are looking for methylation enriched regions. (default: False)')
    argparser.add_argument('--output_label', type=str, default='subCDR', help='Label to use for name column of priorCDR BED file. (default: "subCDR")')

    argparser.add_argument('--save_intermediates', action='store_true', default=False, help="Set to true if you would like to save intermediates(filtered beds+window means). (default: False)")
    
    args = argparser.parse_args()

    sat_types = [st.strip() for st in args.sat_type.split(',')]

    parseCDRBeds = bed_parser(
            mod_code=args.mod_code,
            sat_type=sat_types,
            bedgraph=args.bedgraph,
            min_valid_cov=args.min_valid_cov
        )

    methylation_chrom_dict = parseCDRBeds.process_files(
        bedmethyl_path=args.bedmethyl_path, 
        censat_path=args.censat_path
    )

    priors = calculate_matrices(
        window_size=args.window_size,
        step_size=args.window_size,
        min_prior_size=args.min_prior_size,
        enrichment=args.enrichment,
        output_label=args.output_label
    )

    cdr_priors, window_means = priors.priors_all_chromosomes(methylation_chrom_dict, percentile=args.prior_use_percentile, threshold=args.prior_threshold)


if __name__ == "__main__":
    main()