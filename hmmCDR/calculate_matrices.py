import pandas as pd
import numpy as np 
import pybedtools
import argparse
import os
import concurrent.futures

from hmmCDR.bed_parser import bed_parser

class calculate_matrices:
    def __init__(self, window_size, step_size, min_prior_size,
                 percentile_emissions, enrichment, output_label,
                 w, x, y, z):

        self.window_size = window_size
        self.step_size = step_size
        self.min_prior_size = min_prior_size

        self.percentile_emissions = percentile_emissions

        self.enrichment = enrichment
        self.output_label = output_label

        self.w, self.x, self.y, self.z = w, x, y, z 
        
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
    
    def assign_priors(self, methylation_bedtool, cdr_prior_bedtool):
        """
        Create a DataFrame of genomic regions filtered by mean values and merged based on conditions.

        Args:
            methylation_bedtool (pd.pybedtools.BedTool): BedTool object containing CpG positions and their methylation.
            cdr_priors (pd.pybedtools.BedTool): BedTool object containing regions to be used as priors.
            
        Returns:
            pd.Dataframe: Pandas dataframe chrom positions, fraction modified, and prior state.
        """

        methylation_w_priors_bedtool = methylation_bedtool.intersect(
            cdr_prior_bedtool, 
            c=True
        )

        methylation_w_priors_df = methylation_w_priors_bedtool.saveas().to_dataframe()
        methylation_w_priors_df.columns = ['chrom', 'start', 'end', 'fractionmodified', 'prior']

        return methylation_w_priors_df
    
    def calculate_emission_thresholds(self, methylation_w_priors_df):
        """
        Calculate emission thresholds more efficiently.
        
        Args:
            methylation_w_priors_df (pd.DataFrame): DataFrame with methylation values
            
        Returns:
            list: Sorted emission thresholds
        """
        if not self.percentile_emissions:
            return sorted([self.w, self.x, self.y, self.z])
        
        # Vectorized conversion and filtering in one step
        methylation_scores = pd.to_numeric(
            methylation_w_priors_df['name'], 
            errors='coerce'
        )
        
        # Filter non-zero, non-NaN values and prepend 0
        valid_scores = methylation_scores[methylation_scores > 0]
        scores_with_zero = np.concatenate(([0], valid_scores))
        
        return np.percentile(scores_with_zero, q=[self.w, self.x, self.y, self.z]).tolist()

    def assign_emissions(self, methylation_w_priors_df, emission_thresholds):
        """
        Assign emission states using vectorized operations.
        
        Args:
            bed4Methyl (pd.DataFrame): DataFrame with methylation values
            emission_thresholds (list): Sorted list of threshold values
            
        Returns:
            pd.DataFrame: DataFrame with emissions added
        """
        # Convert values to numeric array
        values = methylation_w_priors_df.iloc[:, 3].values
        
        # Create a matrix of boolean comparisons
        comparison_matrix = values[:, np.newaxis] <= np.array(emission_thresholds)
        
        # Get the first True index for each row (vectorized operation)
        methylation_w_emission_priors_df = methylation_w_priors_df.copy()
        methylation_w_emission_priors_df['emission'] = np.argmax(comparison_matrix, axis=1)
        
        # Handle cases where all comparisons are False
        max_emission = len(emission_thresholds)
        methylation_w_emission_priors_df.loc[~comparison_matrix.any(axis=1), 'emission'] = max_emission - 1

        return methylation_w_emission_priors_df

    def calculate_emission_matrix(self, methylation_w_emission_priors_df):
        """
        Calculate emission matrix using vectorized NumPy operations.

        Computes the probability distribution of emissions for each prior state.

        Args:
            methylation_w_emission_priors_df (pd.DataFrame): DataFrame with prior and emissions columns

        Returns:
            np.ndarray: 2x4 emission matrix with probabilities of emissions for each prior state
        """
        # Create contingency table using NumPy for faster computation
        emission_matrix = np.zeros((2, 4))
        
        for state in [0, 1]:
            state_subset = methylation_w_emission_priors_df[methylation_w_emission_priors_df['prior'] == state]
            emission_counts = np.bincount(state_subset['emission'], minlength=4)
            
            # Normalize to get probabilities
            emission_matrix[state, :] = emission_counts / emission_counts.sum() if emission_counts.sum() > 0 else 0

        # Add a check to ensure rows sum to one (with a small tolerance for floating-point precision)
        assert np.allclose(emission_matrix.sum(axis=1), 1.0, rtol=1e-5), \
            f"Transition matrix rows do not sum to one:\n{emission_matrix}"

        return emission_matrix

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

        # Add a check to ensure rows sum to one (with a small tolerance for floating-point precision)
        assert np.allclose(transition_matrix.sum(axis=1), 1.0, rtol=1e-5), \
            f"Transition matrix rows do not sum to one:\n{transition_matrix}"

        return transition_matrix

    def priors_single_chromosome(self, chrom, methylation_per_chrom, prior_percentile, prior_threshold):
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
        methylation_bedtool = pybedtools.BedTool.from_dataframe(methylation_per_chrom) # Convert the input dataFrame to a bedTool object
        window_bedtool = self.create_chrom_windows(methylation_bedtool) # Create windows
        window_mean_bedtool = self.mean_within_windows(methylation_bedtool, window_bedtool) # Find mean within windows

        # Determine current CDR threshold
        if prior_percentile: 
            prior_threshold = self.find_prior_percentile(window_mean_bedtool, prior_threshold)
        else:
            prior_threshold = prior_threshold

        # Find priors using current threshold
        cdr_prior_bedtool = self.find_priors(window_mean_bedtool, prior_threshold)

        cdr_prior_df = cdr_prior_bedtool.saveas().to_dataframe()
 
        if cdr_prior_df.empty:
            # If no priors found, adjust threshold
            print(f'No prior CDRs Found on {chrom} with prior threshold - {prior_threshold}')
            print(f'Continuing with defaults:\n emission matrix: [[0.002,0.10,0.28,0.60],[0.05,0.85,0.08,0.02]]\ntransition matrix: [[0.9999,0.003],[0.0001,0.997]]')
            emission_matrix = np.array([[0.002,0.10,0.28,0.60],[0.05,0.85,0.08,0.02]])
            transition_matrix = np.array([[0.9999,0.003],[0.0001,0.997]])

            return chrom, cdr_prior_df, cdr_prior_bedtool, emission_matrix, transition_matrix

        # add priors on to methylation bedtool
        methylation_w_priors_df = self.assign_priors(methylation_bedtool, pybedtools.BedTool.from_dataframe(cdr_prior_df))

        # assign emissions to methylation bedtool
        methylation_w_emission_priors_df = self.assign_emissions(methylation_w_priors_df, self.calculate_emission_thresholds(methylation_w_priors_df))

        methylation_w_emission_priors_df.to_csv('methylation_w_emission_priors.tsv',sep='\t',index=None)

        # calculate emission and transition matrices with assigned priors
        emission_matrix = self.calculate_emission_matrix(methylation_w_emission_priors_df)
        transition_matrix = self.calculate_transition_matrix(methylation_w_priors_df)

        return chrom, cdr_prior_df, window_mean_bedtool.saveas().to_dataframe(), methylation_w_emission_priors_df, emission_matrix, transition_matrix

    def priors_all_chromosomes(self, methylation_chrom_dict, prior_percentile, prior_threshold):
        priors_chrom_dict = {}
        windowmean_chrom_dict = {}
        labelled_methylation_chrom_dict = {}
        emission_matrix_chrom_dict = {}
        transition_matrix_chrom_dict = {}
        chromosomes = methylation_chrom_dict.keys()

        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = {
                executor.submit(
                    self.priors_single_chromosome, 
                    chrom,
                    methylation_chrom_dict[chrom], 
                    prior_percentile, 
                    prior_threshold
                ): chrom for chrom in chromosomes
            }

            for future in concurrent.futures.as_completed(futures):
                chrom, cdr_prior, window_mean, labelled_methylation, emission_matrix, transition_matrix = future.result()
                priors_chrom_dict[chrom] = cdr_prior
                windowmean_chrom_dict[chrom] = window_mean
                labelled_methylation_chrom_dict[chrom] = labelled_methylation
                emission_matrix_chrom_dict[chrom] = emission_matrix
                transition_matrix_chrom_dict[chrom] = transition_matrix

        return priors_chrom_dict, windowmean_chrom_dict, labelled_methylation_chrom_dict, emission_matrix_chrom_dict, transition_matrix_chrom_dict


def main():
    argparser = argparse.ArgumentParser(description='Process bedMethyl and CenSat BED file to produce hmmCDR priors')
    
    # required inputs
    argparser.add_argument('bedmethyl', type=str, help='Path to the bedmethyl file')
    argparser.add_argument('censat', type=str, help='Path to the censat BED file')
    argparser.add_argument('output', type=str, help='Path to the output BED file')
    
    # bed_parser arguments
    argparser.add_argument('-m', '--mod_code', type=str, default='m', help='Modification code to filter bedMethyl file (default: "m")')
    argparser.add_argument('-s', '--sat_type', type=str, default='active_hor', help='Comma-separated list of satellite types/names to filter CenSat bed file. (default: "H1L")')
    argparser.add_argument('--bedgraph', action='store_true', default=False, help='Flag indicating if the input is a bedgraph. (default: False)')
    argparser.add_argument('--min_valid_cov', type=int, default=10, help='Minimum valid coverage to consider a methylation site(read from full modkit pileup files). (default: 10)')
    
    # calculate_matrices arguments
    argparser.add_argument('--window_size', type=int, default=1190, help='Window size to calculate prior regions. (default: 1190)')
    argparser.add_argument('--step_size', type=int, default=1190, help='Step size when calculation windows for priors. (default: 1190)')
    argparser.add_argument('--prior_threshold', type=float, default=30.0, help='Threshold for determining if a window is a CDR. Uses this percentile if --prior_use_percentile is passed (default: 30.0)')
    argparser.add_argument('--prior_use_percentile', action='store_true', default=False, help='Whether or not to use percentile when calculating windowing priors. (default: False)')
    argparser.add_argument('--min_prior_size', type=int, default=8330, help='Minimum size for CDR regions. (default: 8330)')
    argparser.add_argument('--enrichment', action='store_true', default=False, help='Enrichment flag. Pass in if you are looking for methylation enriched regions. (default: False)')
    argparser.add_argument('--output_label', type=str, default='subCDR', help='Label to use for name column of priorCDR BED file. (default: "subCDR")')
    argparser.add_argument('--percentile_emissions', action='store_true', default=False, help='Use values for flags w,x,y,z as raw threshold cutoffs for each emission category. (default: False)')
    argparser.add_argument('-w', type=float, default=0.0, help='Threshold of non-zero methylation percentile to be classified as None (default: 0.0)')
    argparser.add_argument('-x', type=float, default=33.3, help='Threshold of non-zero methylation percentile to be classified as low (default: 33.3)')
    argparser.add_argument('-y', type=float, default=66.6, help='Threshold of non-zero methylation percentile to be classified as medium (default: 66.6)')
    argparser.add_argument('-z', type=float, default=100.0, help='Threshold of non-zero methylation percentile to be classified as high (default: 100.0)')

    args = argparser.parse_args()
    sat_types = [st.strip() for st in args.sat_type.split(',')]
    output_prefix = os.path.splitext(args.output)[0]

    parseCDRBeds = bed_parser(
        mod_code = args.mod_code,
        sat_type = sat_types,
        bedgraph = args.bedgraph,
        min_valid_cov = args.min_valid_cov
    )

    methylation_chrom_dict = parseCDRBeds.process_files(
        bedmethyl_path = args.bedmethyl, 
        censat_path = args.censat
    )

    priors = calculate_matrices(
        window_size = args.window_size,
        step_size = args.window_size,
        min_prior_size = args.min_prior_size,
        enrichment = args.enrichment,
        percentile_emissions = args.percentile_emissions,
        w = args.w, x = args.x, y = args.y, z = args.z,
        output_label = args.output_label
    )

    priors_chrom_dict, windowmean_chrom_dict, labelled_methylation_chrom_dict, emission_matrix_chrom_dict, transition_matrix_chrom_dict = priors.priors_all_chromosomes(
        methylation_chrom_dict = methylation_chrom_dict, 
        prior_percentile = args.prior_use_percentile, 
        prior_threshold = args.prior_threshold
    )

    # Concatenate, sort and save bed files
    combined_priors = pd.concat(priors_chrom_dict.values(), ignore_index=True)
    combined_priors = combined_priors.sort_values(by=combined_priors.columns[:2].tolist())
    combined_priors.to_csv(f'{output_prefix}_priors.bed', sep='\t', index=False, header=False)

    combined_windowmean = pd.concat(windowmean_chrom_dict.values(), ignore_index=True)
    combined_windowmean = combined_windowmean.sort_values(by=combined_windowmean.columns[:2].tolist())
    combined_windowmean.to_csv(f'{output_prefix}_windowmean.bedgraph', sep='\t', index=False, header=False)

    # Save matrices to a text file with key: numpy array representation
    with open(f'{output_prefix}_emission_matrices.txt', 'w') as f:
        for key, matrix in emission_matrix_chrom_dict.items():
            f.write(f"{key}: {matrix.tolist()}\n")

    with open(f'{output_prefix}_transition_matrices.txt', 'w') as f:
        for key, matrix in transition_matrix_chrom_dict.items():
            f.write(f"{key}: {matrix.tolist()}\n")


if __name__ == "__main__":
    main()