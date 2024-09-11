import pandas as pd
import numpy as np 
import pybedtools
import argparse
import os
import concurrent.futures

from hmmCDR.parser import hmmCDR_parser

class hmmCDR_prior_finder:
    def __init__(self, merge_distance, window_size, min_size, prior_percentile,  enrichment, output_label):
        self.merge_distance = merge_distance
        self.window_size = window_size
        self.min_size = min_size

        self.prior_percentile = prior_percentile

        self.enrichment = enrichment
        self.output_label = output_label
        
        self.retries = 0
        

    def create_windows(self, bed4Methyl):
        chrom = bed4Methyl.iloc[0]['chrom']
        min_val, max_val = bed4Methyl['start'].min(), bed4Methyl['end'].max()
        windows = pd.DataFrame(
            [[chrom, start, start + self.window_size] 
            for start in range(int(min_val), int(max_val) - self.window_size + 1, self.window_size)],
            columns=[0, 1, 2]
        )
        return windows

    def mean_within_windows(self, bed4Methyl, windows):
        return pybedtools.BedTool.from_dataframe(windows).map(pybedtools.BedTool.from_dataframe(bed4Methyl), c=4, o='mean').to_dataframe(names=['chrom', 'start', 'end', 'mean_value'])

    def calculate_prior_percentiles(self, window_means):
        mean_values = pd.to_numeric(window_means['mean_value'].replace('.', np.nan)).dropna()
        return ( np.percentile(mean_values, self.prior_percentile) )

    def create_prior_dataframe(self, merge_distance, windows_mean_df, priorCDR_score):
        windows_mean_df['mean_value'] = pd.to_numeric(windows_mean_df['mean_value'], errors='coerce')
        windows_mean_df = windows_mean_df.dropna(subset=['mean_value'])
        if self.enrichment:
            condition = windows_mean_df['mean_value'] > priorCDR_score
        else:
            condition = windows_mean_df['mean_value'] < priorCDR_score

        priorCDR_windows = windows_mean_df[condition]

        merged_windows = pybedtools.BedTool.from_dataframe(priorCDR_windows).merge(d=merge_distance)
        merged_df = merged_windows.to_dataframe()

        if merged_df.empty:
            return pd.DataFrame()
        
        return merged_df[merged_df['end'] - merged_df['start'] >= self.min_size]

    def combine_beds(self, priorCDR_df, priorTransition_df):
        priorCDR_df[3], priorTransition_df[3] = self.output_label, f"{self.output_label}_transition"
        return pd.concat([priorCDR_df, priorTransition_df], ignore_index=True).sort_values(by=1).reset_index(drop=True)

    def priors_single_chromosome(self, chrom, bed4Methyl_chrom):
        windows = self.create_windows(bed4Methyl_chrom)
        windows_mean = self.mean_within_windows(bed4Methyl_chrom, windows)
        cdr_threshold = self.calculate_prior_percentiles(windows_mean)
        prior_cdrs = self.create_prior_dataframe(self.merge_distance, windows_mean, cdr_threshold)
        
        if prior_cdrs.empty:
            print(f'No prior subCDRs Detected for {chrom} with percentile - {self.prior_percentile}')
            
            if self.retries < 25:
                self.prior_percentile = min(self.prior_percentile + 1, 100)

                self.retries += 1
                print(f'Retrying with Prior Percentile = {self.prior_percentile}')
                print(f'This is retry attempt {self.retries}/25.')

                return self.priors_single_chromosome(chrom, bed4Methyl_chrom)
            raise RuntimeError(f"Failed to detect prior subCDRs for {chrom} after {self.retries} retries with final percentile - {self.prior_percentile}")

        return chrom, prior_cdrs

    def priors_all_chromosomes(self, bed4Methyl_chrom_dict):
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
    argparser.add_argument('-s', '--sat_type', type=str, default='H1L', help='Comma-separated list of satellite types/names to filter CenSat bed file. (default: "H1L")')
    argparser.add_argument('--bedgraph', action='store_true', help='Flag indicating if the input is a bedgraph. (default: False)')
    argparser.add_argument('--min_valid_cov', type=int, default=10, help='Minimum valid coverage to consider a methylation site(read from full modkit pileup files). (default: 10)')
    
    # hmmCDRprior Arguments
    argparser.add_argument('--window_size', type=int, default=1020, help='Window size to calculate prior regions. (default: 1020)')
    argparser.add_argument('--prior_percentile', type=int, default=5, help='Percentile for finding CDR prior regions. (default: 5)')
    argparser.add_argument('--merge_distance', type=int, default=1021, help='Distance to merge adjacent CDR regions. (default: 1021)')
    argparser.add_argument('--min_size', type=int, default=3000, help='Minimum size for CDR regions. (default: 3000)')
    argparser.add_argument('--enrichment', action='store_true', help='Enrichment flag. Pass in if you are looking for methylation enriched regions. (default: False)')
    argparser.add_argument('--output_label', type=str, default='subCDR', help='Label to use for name column of priorCDR BED file. (default: "subCDR")')

    argparser.add_argument('--save_intermediates', action='store_true', default=False, help="Set to true if you would like to save intermediates(filtered beds+window means). (default: False)")
    
    args = argparser.parse_args()

    sat_types = [st.strip() for st in args.sat_type.split(',')]

    CDRparser = hmmCDRparse(
            mod_code=args.mod_code,
            sat_type=sat_types,
            bedgraph=args.bedgraph,
            min_valid_cov=args.min_valid_cov
        )

    bed4Methyl_chrom_dict, cenSat_chrom_dict = CDRparser.process_files(
        bedMethyl_path=args.bedMethyl_path, 
        cenSat_path=args.cenSat_path
    )

    priors = hmmCDRprior(
        window_size=args.window_size,
        min_size=args.min_size,
        prior_percentile=args.prior_percentile,
        merge_distance=args.merge_distance, 
        enrichment=args.enrichment,
        output_label=args.output_label
    )

    concatenated_priors = pd.concat(
        priors.priors_all_chromosomes(bed4Methyl_chrom_dict).values(), axis=0
    )
    
    concatenated_priors.to_csv(args.output_path, sep='\t', index=False, header=False)
    print(f"hmmCDR Priors saved to: {args.output_path}")


if __name__ == "__main__":
    main()