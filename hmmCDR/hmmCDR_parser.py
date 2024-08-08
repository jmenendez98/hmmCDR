import pandas as pd
import numpy as np 
import pybedtools
import argparse
import os


class hmmCDR_parser:
    '''
    CLASS DOCSTRING
    '''
    def __init__(self, bedMethyl_path, cenSat_path, 
                 mod_code, sat_type, bedgraph=False):
        '''
        INIT DOCSTRING
        '''
        self.mod_code = mod_code
        self.sat_type = sat_type
        self.bedgraph = bedgraph

        bedMethyl_df = pd.read_csv(bedMethyl_path, sep='\t', header=None)
        cenSat_df = pd.read_csv(cenSat_path, sep='\t', header=None)

        if not bedgraph:
            if len(bedMethyl_df.columns) != 18:
                raise ValueError("Valid bedMethyl must have 18 columns when the bedgraph flag is not passed.")
            self.bedgraphMethyl = self.filter_bedMethyl(bedMethyl_df, self.mod_code)
        else:
            if len(bedMethyl_df.columns) != 4:
                raise ValueError("Valid bedgraphMethyl must have 4 columns when the bedgraph flag is passed.")
            self.bedgraphMethyl = bedMethyl_df

        self.subset_cenSat = self.filter_cenSat(cenSat_df, self.sat_type)

        self.subset_bedgraphMethyl = self.intersect_files(self.bedgraphMethyl, self.subset_cenSat)

    def filter_bedMethyl(self, bedMethyl, mod_code):
        '''
        DOCSTRING
        '''
        filtered_bedMethyl = bedMethyl[bedMethyl[3] == mod_code]
        bedgraphMethyl = filtered_bedMethyl.iloc[:, [0, 1, 2, 10]]
        return bedgraphMethyl

    def filter_cenSat(self, cenSat, sat_type):
        '''
        DOCSTRING
        '''
        filtered_cenSat = cenSat[cenSat[3].str.contains(sat_type)]
        return filtered_cenSat

    def intersect_files(self, bedgraphMethyl, subset_cenSat):
        '''
        DOCSTRING
        '''
        bedMethyl_bedtool = pybedtools.BedTool.from_dataframe(bedgraphMethyl)
        regions_bedtool = pybedtools.BedTool.from_dataframe(subset_cenSat)
        intersected = bedMethyl_bedtool.intersect(regions_bedtool, wa=True, u=True)
        intersected_df = intersected.to_dataframe(names=[0, 1, 2, 3])
        return intersected_df
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process bedMethyl and CenSat BED file to produce hmmCDR priors')
    # Required arguments
    parser.add_argument('bedMethyl_path', type=str, help='Path to the bedMethyl file')
    parser.add_argument('cenSat_path', type=str, help='Path to the CenSat BED file')
    parser.add_argument('output_prefix', type=str, help='Path to the output priorCDRs BED file')
    # Optional arguments with default values
    parser.add_argument('--bedgraph', action='store_true', help='Flag indicating if the input is a bedgraph. (default: False)')
    parser.add_argument('-m', '--mod_code', type=str, default='m', help='Modification code to filter bedMethyl file (default: "m")')
    parser.add_argument('-s', '--sat_type', type=str, default='H1L', help='Satellite type/name to filter CenSat bed file. (default: "H1L")')
    args = parser.parse_args()

    hmmCDRparser = hmmCDR_parser(bedMethyl_path=args.bedMethyl_path,
                           cenSat_path=args.cenSat_path,
                           mod_code=args.mod_code,
                           sat_type=args.sat_type,
                           bedgraph=args.bedgraph)
    
    # Save the filtered DataFrames to TSV files
    hmmCDRparser.subset_bedgraphMethyl.to_csv(f'{args.output_prefix}_filtered_bedMethyl.bedgraph', sep='\t', index=False, header=False)
    hmmCDRparser.subset_cenSat.to_csv(f'{args.output_prefix}_filtered_cenSat.bed', sep='\t', index=False, header=False)
    print(f"Filtered bedgraphMethyl saved to: {args.output_prefix}_filtered_bedgraphMethyl.bedgraph")
    print(f"Filtered CenSat saved to: {args.output_prefix}_filtered_cenSat.bed")