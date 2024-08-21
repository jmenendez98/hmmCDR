import pandas as pd
import numpy as np 
import pybedtools
import argparse
import os
import concurrent.futures

from hmmCDR.hmmCDRparse import hmmCDRparse

class hmmCDRprior:
    '''
    The hmmCDR_priors class processes bedMethyl data to identify candidate
    CDR (Centromere Dig Regions) and transition regions. It calculates
    these regions based on user-defined parameters such as window size,
    percentiles, and minimum region size. The results are combined into
    a BED file format.
    '''
    def __init__(self, window_size=1020, minCDR_size=3000,
                 priorCDR_percent=5, priorTransition_percent=10, 
                 enrichment=False, output_label='CDR'):
        '''
        Initializes the hmmCDR_priors class with the given parameters and
        processes the data to identify CDR and transition regions.

        Parameters:
        - bedgraphMethyl (pd.DataFrame): DataFrame containing the bedMethyl data.
        - output_path (str): Path where the output BED file will be saved.
        - window_size (int): Size of the windows used to calculate prior regions. Default is 1020.
        - minCDR_size (int): Minimum size for CDR regions. Default is 3000.
        - priorCDR_percent (int): Percentile for identifying prior CDR regions. Default is 5.
        - priorTransition_percent (int): Percentile for identifying prior transition regions. Default is 10.
        - output_label (str): Label used in the output BED file for CDR regions. Default is 'CDR'.
        - enrichment (bool): Flag to indicate if looking for methylation-enriched regions. Default is False.
        - save_intermediates (bool): Flag to indicate if intermediate files should be saved. Default is False.
        '''

        # all hmmCDR_prior parameters can be optionally changed
        self.window_size = window_size
        self.minCDR_size = minCDR_size
        self.priorCDR_percent = priorCDR_percent
        self.priorTransition_percent = priorTransition_percent
        self.enrichment = enrichment
        self.output_label = output_label

        self.retries = 0
        

    def create_windows(self, bed4Methyl):
        '''
        Creates a DataFrame of genomic windows based on the bedMethyl data.

        Parameters:
        - bed4Methyl (pd.DataFrame): DataFrame containing intersected bedMethyl data.

        Returns:
        - pd.DataFrame: DataFrame containing the created windows.
        '''
        min_val = int(bed4Methyl['start'].min())
        max_val = int(bed4Methyl['end'].max())
        regions = []
        current_start = min_val
        while current_start + self.window_size <= max_val:
            current_end = current_start + self.window_size
            chrom = bed4Methyl.iloc[0]['chrom'] 
            regions.append([chrom, current_start, current_end])
            current_start = current_end
        windows = pd.DataFrame(regions, columns=[0, 1, 2])
        return windows

    def mean_within_windows(self, bed4Methyl, windows):
        '''
        Calculates the mean methylation within each window.

        Parameters:
        - bedMethyl_df (pd.DataFrame): DataFrame containing the bedMethyl data.
        - windows_df (pd.DataFrame): DataFrame containing the genomic windows.

        Returns:
        - pd.DataFrame: DataFrame containing the mean methylation values within each window.
        '''
        bed4Methyl_bedtool = pybedtools.BedTool.from_dataframe(bed4Methyl)
        windows_bedtool = pybedtools.BedTool.from_dataframe(windows)
        window_means_map = windows_bedtool.map(bed4Methyl_bedtool, c=4, o='mean')
        window_means = window_means_map.to_dataframe(names=[0, 1, 2, 'mean_value'])
        return window_means

    def calculate_prior_percentiles(self, window_means):
        '''
        Calculates the percentile scores for CDR and transition regions.

        Parameters:
        - windows_mean_df (pd.DataFrame): DataFrame containing the mean methylation values within each window.

        Returns:
        - tuple: A tuple containing the CDR score and transition score.
        '''
        window_means['mean_value'].replace('.', np.nan, inplace=True)
        mean_values = window_means['mean_value'].dropna().astype(float)
        priorCDR_score = np.percentile(mean_values, q=self.priorCDR_percent)
        priorTransition_score = np.percentile(mean_values, q=self.priorTransition_percent)
        return priorCDR_score, priorTransition_score

    def create_priorCDR_dataframe(self, windows_mean_df, priorCDR_score):
        '''
        Creates a DataFrame of prior CDR regions based on the calculated score.

        Parameters:
        - windows_mean_df (pd.DataFrame): DataFrame containing the mean methylation values within each window.
        - priorCDR_score (float): Score used to identify CDR regions.

        Returns:
        - pd.DataFrame: DataFrame containing the prior CDR regions.
        '''
        windows_mean_df['mean_value'] = pd.to_numeric(windows_mean_df['mean_value'], errors='coerce')
        windows_mean_df = windows_mean_df.dropna(subset=['mean_value'])
        if self.enrichment:
            windows_below_priorCDR_score = windows_mean_df[windows_mean_df['mean_value'] > priorCDR_score]
        else:
            windows_below_priorCDR_score = windows_mean_df[windows_mean_df['mean_value'] < priorCDR_score]
        windows_bedtool = pybedtools.BedTool.from_dataframe(windows_below_priorCDR_score)
        merged_windows = windows_bedtool.merge()
        merged_df = merged_windows.to_dataframe(names=[0, 1, 2])
        merged_df['size'] = merged_df[2] - merged_df[1]
        filtered_merged_df = merged_df[merged_df['size'] >= self.minCDR_size]
        filtered_merged_df = filtered_merged_df.drop(columns=['size'])
        return filtered_merged_df

    def create_priorTransition_dataframe(self, windows_mean_df, priorTransition_score, CDR_df):
        '''
        Creates a DataFrame of prior transition regions based on the calculated score and CDR regions.

        Parameters:
        - windows_mean_df (pd.DataFrame): DataFrame containing the mean methylation values within each window.
        - priorTransition_score (float): Score used to identify transition regions.
        - CDR_df (pd.DataFrame): DataFrame containing the prior CDR regions.

        Returns:
        - pd.DataFrame: DataFrame containing the prior transition regions.
        '''
        windows_mean_df['mean_value'] = pd.to_numeric(windows_mean_df['mean_value'], errors='coerce')
        windows_mean_df = windows_mean_df.dropna(subset=['mean_value'])
        if self.enrichment:
            windows_below_priorTransition_score = windows_mean_df[windows_mean_df['mean_value'] > priorTransition_score]
        else:
            windows_below_priorTransition_score = windows_mean_df[windows_mean_df['mean_value'] < priorTransition_score]
        windows_bedtool = pybedtools.BedTool.from_dataframe(windows_below_priorTransition_score)
        merged_windows = windows_bedtool.merge()
        CDR_bedtool = pybedtools.BedTool.from_dataframe(CDR_df)
        intersected_bedtool = merged_windows.intersect(CDR_bedtool, wa=True, wb=True)
        intersected_df = intersected_bedtool.to_dataframe(names=[0, 1, 2, 'CDR_1', 'CDR_2', 'CDR_3', 'CDR_4'])
        intersected_df = intersected_df[[0, 1, 2]]
        intersected_bedtool = pybedtools.BedTool.from_dataframe(intersected_df)
        priorTransition_bedtool = intersected_bedtool.subtract(CDR_bedtool)
        priorTransition_df = priorTransition_bedtool.to_dataframe(names=[0, 1, 2])
        return priorTransition_df
    
    def combine_beds(self, priorCDR_df, priorTransition_df):
        '''
        Combines the CDR and transition regions into a single DataFrame.

        Parameters:
        - priorCDR_df (pd.DataFrame): DataFrame containing the prior CDR regions.
        - priorTransition_df (pd.DataFrame): DataFrame containing the prior transition regions.

        Returns:
        - pd.DataFrame: Combined DataFrame with both CDR and transition regions.
        '''
        priorCDR_df[3] = self.output_label
        priorTransition_df[3] = self.output_label + "_transition"
        combined_df = pd.concat([priorCDR_df, priorTransition_df], ignore_index=True)
        combined_df = combined_df.sort_values(by=1).reset_index(drop=True)
        return combined_df

    def priors_single_chromosome(self, chrom, bed4Methyl_chrom):
        # Prior processing specific to each chromosome
        windows = self.create_windows(bed4Methyl_chrom)
        windows_mean = self.mean_within_windows(bed4Methyl_chrom, windows)
        cdr_score, transition_score = self.calculate_prior_percentiles(windows_mean)
        priorCDRs = self.create_priorCDR_dataframe(windows_mean, cdr_score)
        priorTranstions = self.create_priorTransition_dataframe(windows_mean, transition_score, priorCDRs)
        hmmCDRpriors = self.combine_beds(priorCDRs, priorTranstions)
        if hmmCDRpriors.empty:
            print(f'No Priors Detected for {chrom} with settings: CDR Percentile - {self.priorCDR_percent}, Transition Percentile - {self.priorTransition_percent}')
            
            if self.retries < 25:
                if 0 < self.priorCDR_percent < 100:
                    self.priorCDR_percent += 1
                if 0 < self.priorTransition_percent < 100:
                    self.priorTransition_percent += 1
                self.retries += 1
                print(f'Retrying with CDR Percentile = {self.priorCDR_percent}, Transition Percentile = {self.priorTransition_percent}')
                print(f'This is the {self.retries}/25 retry attempt.')
                return self.priors_single_chromosome(chrom, bed4Methyl_chrom)
            else:
                raise RuntimeError(f"Failed to detect priors for {chrom} after {self.retries} retries with final settings: CDR Percentile - {self.priorCDR_percent}, Transition Percentile - {self.priorTransition_percent}")

        return chrom, hmmCDRpriors

    def priors_all_chromosomes(self, bed4Methyl_chrom_dict):
        '''
        Processes all chromosomes in parallel using concurrent futures.

        Parameters:
        -----------
        bed4Methyl_chrom_dict : dict
            A dictionary with chromosome names as keys and DataFrames (from bedMethyl) as values.
        
        cenSat : pd.DataFrame
            The DataFrame containing cenSat annotations.

        Returns:
        --------
        dict
            Two dictionaries:
            - final_bed4Methyl_chrom_dict: A dictionary with chromosome names as keys and processed DataFrames as values.
            - cenSat_chrom_dict: A dictionary with chromosome names as keys and filtered cenSat DataFrames as values.
        '''
        hmmCDRpriors_chrom_dict = {}
        chromosomes = bed4Methyl_chrom_dict.keys()

        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = {
                executor.submit(
                    self.priors_single_chromosome, chrom,
                    bed4Methyl_chrom_dict[chrom]
                ): chrom for chrom in chromosomes
            }

            for future in concurrent.futures.as_completed(futures):
                chrom, priors_df = future.result()

                hmmCDRpriors_chrom_dict[chrom] = priors_df

        self.chromosomes = chromosomes
        self.hmmCDRpriors_chrom_dict = hmmCDRpriors_chrom_dict

        return hmmCDRpriors_chrom_dict


def main():
    argparser = argparse.ArgumentParser(description='Process bedMethyl and CenSat BED file to produce hmmCDR priors')
    
    # Required arguments
    argparser.add_argument('bedMethyl_path', type=str, help='Path to the bedMethyl file')
    argparser.add_argument('cenSat_path', type=str, help='Path to the CenSat BED file')
    argparser.add_argument('output_path', type=str, help='Path to the output priorCDRs BED file')
    
    # Parser Arguments
    argparser.add_argument('-m', '--mod_code', type=str, default='m', help='Modification code to filter bedMethyl file (default: "m")')
    argparser.add_argument('-s', '--sat_type', type=str, default='H1L', help='Satellite type/name to filter CenSat bed file. (default: "H1L")')
    argparser.add_argument('--bedgraph', action='store_true', help='Flag indicating if the input is a bedgraph. (default: False)')
    argparser.add_argument('--rolling_window', type=int, default=0, help='Flag indicating whether or not to use a rolling average and the rolling avg window size. If set to 0 no rolling averages are used. (defualt: 0)')
    argparser.add_argument('--min_valid_cov', type=int, default=10, help='Minimum Valid Coverage to consider a methylation site. (default: 10)')
    
    # hmmCDRprior Arguments
    argparser.add_argument('-w', '--window_size', type=int, default=1020, help='Window size to calculate prior regions. (default: 1020)')
    argparser.add_argument('--priorCDR_percent', type=int, default=5, help='Percentile for finding priorCDR regions. (default: 5)')
    argparser.add_argument('--priorTransition_percent', type=int, default=10, help='Percentile for finding priorTransition regions. (default: 10)')
    argparser.add_argument('--minCDR_size', type=int, default=3000, help='Minimum size for CDR regions. (default: 3000)')
    argparser.add_argument('--enrichment', action='store_true', help='Enrichment flag. Pass in if you are looking for methylation enriched regions. (default: False)')
    argparser.add_argument('--save_intermediates', action='store_true', default=False, help="Set to true if you would like to save intermediates(filtered beds+window means). (default: False)")
    argparser.add_argument('--output_label', type=str, default='CDR', help='Label to use for name column of priorCDR BED file. (default: "CDR")')
    
    args = argparser.parse_args()

    parser = hmmCDRparse(
        bedMethyl_path=args.bedMethyl_path,
        cenSat_path=args.cenSat_path,
        mod_code=args.mod_code,
        sat_type=args.sat_type,
        rolling_window=args.rolling_window,
        bedgraph=args.bedgraph
    )

    cenSat = parser.read_cenSat(path=parser.cenSat_path)
    bedMethyl = parser.read_bedMethyl(path=parser.bedMethyl_path)
    bed4Methyl_chrom_dict, cenSat_chrom_dict = parser.parse_all_chromosomes(bedMethyl=bedMethyl, cenSat=cenSat)

    priors = hmmCDRprior(
        window_size=args.window_size,
        minCDR_size=args.minCDR_size,
        priorCDR_percent=args.priorCDR_percent,
        priorTransition_percent=args.priorTransition_percent,
        enrichment=args.enrichment,
        output_label=args.output_label
    )

    hmmCDRpriors = priors.priors_all_chromosomes(bed4Methyl_chrom_dict)

    # Combine all chromosomes and save the output
    concatenated_priors = pd.concat(hmmCDRpriors.values(), axis=0)
    concatenated_priors.to_csv(args.output_path, sep='\t', index=False, header=False)
    print(f"hmmCDR Priors saved to: {args.output_path}")


if __name__ == "__main__":
    main()