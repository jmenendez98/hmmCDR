import pandas as pd
import numpy as np 
import pybedtools
import argparse
import os

from hmmCDR_parser import hmmCDR_parser

class hmmCDR_priors:
    def __init__(self, bedgraphMethyl, output_path, 
                 window_size=1020, minCDR_size=3000,
                 priorCDR_percent=5, priorTransition_percent=10, 
                 output_label='CDR', enrichment=False,
                 save_intermediates=False):
        
        # required
        self.bedgraphMethyl = bedgraphMethyl
        self.output_path = output_path

        # optionally changed
        self.window_size = window_size
        self.minCDR_size = minCDR_size
        self.priorCDR_percent = priorCDR_percent
        self.priorTransition_percent = priorTransition_percent
        self.enrichment = enrichment
        self.output_label = output_label
        self.save_intermediates=save_intermediates

        # calculated
        windows = self.create_windows(self.bedgraphMethyl)
        windows_mean = self.mean_within_windows(self.bedgraphMethyl, windows)
        self.cdr_score, self.transition_score = self.calculate_percentiles(windows_mean)

        self.priorCDR = self.create_priorCDR_dataframe(windows_mean, self.cdr_score)
        self.priorTranstion = self.create_priorTransition_dataframe(windows_mean, self.transition_score, self.priorCDR)

        self.hmmCDRpriors = self.combine_beds(self.priorCDR, self.priorTranstion)
        

    def create_windows(self, intersected_bedMethyl):
        min_val = int(intersected_bedMethyl[1].min())
        max_val = int(intersected_bedMethyl[1].max())
        regions = []
        current_start = min_val
        while current_start + self.window_size <= max_val:
            current_end = current_start + self.window_size
            regions.append([intersected_bedMethyl[0][0], current_start, current_end])
            current_start = current_end
        regions_df = pd.DataFrame(regions, columns=[0, 1, 2])
        return regions_df

    def mean_within_windows(self, bedMethyl_df, windows_df):
        intersected_bedtool = pybedtools.BedTool.from_dataframe(bedMethyl_df)
        windows_bedtool = pybedtools.BedTool.from_dataframe(windows_df)
        window_means = windows_bedtool.map(intersected_bedtool, c=4, o='mean')
        window_means_df = window_means.to_dataframe(names=[0, 1, 2, 'mean_value'])
        return window_means_df

    def calculate_percentiles(self, windows_mean_df):
        windows_mean_df['mean_value'].replace('.', np.nan, inplace=True)
        mean_values = windows_mean_df['mean_value'].dropna().astype(float)
        priorCDR_score = np.percentile(mean_values, q=self.priorCDR_percent)
        priorTransition_score = np.percentile(mean_values, q=self.priorTransition_percent)
        return priorCDR_score, priorTransition_score

    def create_priorCDR_dataframe(self, windows_mean_df, priorCDR_score):
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
        priorCDR_df[3] = self.output_label
        priorTransition_df[3] = self.output_label + "_transition"
        combined_df = pd.concat([priorCDR_df, priorTransition_df], ignore_index=True)
        combined_df = combined_df.sort_values(by=1).reset_index(drop=True)
        return combined_df


def main():
    parser = argparse.ArgumentParser(description='Process bedMethyl and CenSat BED file to produce hmmCDR priors')
    
    # Required arguments
    parser.add_argument('bedMethyl_path', type=str, help='Path to the bedMethyl file')
    parser.add_argument('cenSat_path', type=str, help='Path to the CenSat BED file')
    parser.add_argument('output_path', type=str, help='Path to the output priorCDRs BED file')
    
    # Optional arguments with default values
    parser.add_argument('-m', '--mod_code', type=str, default='m', help='Modification code to filter bedMethyl file (default: "m")')
    parser.add_argument('-s', '--sat_type', type=str, default='H1L', help='Satellite type/name to filter CenSat bed file. (default: "H1L")')
    parser.add_argument('--bedgraph', action='store_true', help='Flag indicating if the input is a bedgraph. (default: False)')
    parser.add_argument('-w', '--window_size', type=int, default=1020, help='Window size to calculate prior regions. (default: 1020)')
    parser.add_argument('--priorCDR_percent', type=int, default=5, help='Percentile for finding priorCDR regions. (default: 5)')
    parser.add_argument('--priorTransition_percent', type=int, default=10, help='Percentile for finding priorTransition regions. (default: 10)')
    parser.add_argument('--minCDR_size', type=int, default=3000, help='Minimum size for CDR regions. (default: 3000)')
    parser.add_argument('--enrichment', action='store_true', help='Enrichment flag. Pass in if you are looking for methylation enriched regions. (default: False)')
    parser.add_argument('--save_intermediates', action='store_true', default=False, help="Set to true if you would like to save intermediates(filtered beds+window means). (default: False)")
    parser.add_argument('--output_label', type=str, default='CDR', help='Label to use for name column of priorCDR BED file. (default: "CDR")')
    
    args = parser.parse_args()

    CDRparser = hmmCDR_parser(
        bedMethyl_path=args.bedMethyl_path,
        cenSat_path=args.cenSat_path,
        mod_code=args.mod_code,
        sat_type=args.sat_type,
        bedgraph=args.bedgraph
    )

    CDRpriors = hmmCDR_priors(
        bedgraphMethyl=CDRparser.subset_bedgraphMethyl,
        output_path=args.output_path,
        window_size=args.window_size,
        minCDR_size=args.minCDR_size,
        priorCDR_percent=args.priorCDR_percent,
        priorTransition_percent=args.priorTransition_percent,
        enrichment=args.enrichment,
        save_intermediates=args.save_intermediates,
        output_label=args.output_label
    )

    CDRpriors.hmmCDRpriors.to_csv(args.output_path, sep='\t', header=False, index=False)

if __name__ == "__main__":
    main()