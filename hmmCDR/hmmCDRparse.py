import pandas as pd
import numpy as np 
import pybedtools
import argparse
import os


class hmmCDR_parser:
    '''
    This class processes bedMethyl and CenSat BED files to produce HMM CDR priors. 
    It filters and intersects the input files based on specified modification codes 
    and satellite types. The resulting filtered data is used to create intersected 
    regions, which can be saved as output files.
    
    Attributes:
    -----------
    mod_code : str
        Modification code to filter the bedMethyl file.
    sat_type : str
        Satellite type/name to filter the CenSat BED file.
    bedgraph : bool
        Flag indicating whether the input file is in bedgraph format.
    bedgraphMethyl : DataFrame
        The filtered bedMethyl DataFrame.
    subset_cenSat : DataFrame
        The filtered CenSat DataFrame.
    subset_bedgraphMethyl : DataFrame
        The DataFrame containing intersected regions between the filtered 
        bedMethyl and CenSat files.
    '''
    def __init__(self, bedMethyl_path, cenSat_path, 
                 mod_code, sat_type, bedgraph=False):
        '''
        Initializes the hmmCDR_parser class by reading and processing the bedMethyl 
        and CenSat BED files, filtering them according to the provided modification 
        code and satellite type, and intersecting the filtered regions.
        
        Parameters:
        -----------
        bedMethyl_path : str
            Path to the bedMethyl file.
        cenSat_path : str
            Path to the CenSat BED file.
        mod_code : str
            Modification code to filter the bedMethyl file.
        sat_type : str
            Satellite type/name to filter the CenSat BED file.
        bedgraph : bool, optional
            Flag indicating if the input is a bedgraph file. (default is False)
        
        Raises:
        -------
        ValueError
            If the bedMethyl file does not have the expected number of columns 
            based on the bedgraph flag.
        '''
        self.mod_code = mod_code
        self.sat_type = sat_type
        self.bedgraph = bedgraph
        self.bedMethyl_path = bedMethyl_path
        self.cenSat_path = cenSat_path


    def read_bedMethyl(self, path):
        # Check if bedMethyl file exists and is not empty
        if not os.path.exists(path) or os.stat(path).st_size == 0:
            raise ValueError(f"bedMethyl file {path} does not exist or is empty.")
        self.bedMethyl = pd.read_csv(path, sep='\t', header=None)
        return self.bedMethyl
    
    def read_cenSat(self, path):
        # Check if cenSat file exists and is not empty
        if not os.path.exists(path) or os.stat(path).st_size == 0:
            raise ValueError(f"cenSat file {path} does not exist or is empty.")
        self.cenSat= pd.read_csv(path, sep='\t', header=None)
        return self.cenSat

    def filter_bedMethyl(self, bedMethyl, mod_code):
        '''
        Filters the bedMethyl DataFrame by the specified modification code 
        and subsets it to include only relevant columns.

        Parameters:
        -----------
        bedMethyl : DataFrame
            The DataFrame containing the bedMethyl data.
        mod_code : str
            The modification code to filter the bedMethyl data.

        Returns:
        --------
        DataFrame
            A DataFrame containing the filtered bedMethyl data with only 
            the relevant columns (chromosome, start, end, and methylation level).

        Raises:
        -------
        ValueError
            If the filtering results in an empty DataFrame.
        '''
        if self.bedgraph:
            if len(bedMethyl.columns) != 4:
                raise ValueError("Valid bedgraph must have 4 columns when the bedgraph flag is passed.")
            self.bed4Methyl = bedMethyl
            return self.bed4Methyl
        
        if len(bedMethyl.columns) != 18:
                raise ValueError("Valid bedMethyl must have 18 columns when the bedgraph flag is not passed.")

        filtered_bedMethyl = bedMethyl[bedMethyl[3] == mod_code]
        self.bed4Methyl = filtered_bedMethyl.iloc[:, [0, 1, 2, 10]]
        if self.bed4Methyl.empty:
            raise ValueError("Filtering bedMethyl by the specified modification code resulted in an empty DataFrame.")
        return self.bed4Methyl

    def filter_cenSat(self, cenSat, sat_type):
        '''
        Filters the CenSat DataFrame by the specified satellite type.

        Parameters:
        -----------
        cenSat : DataFrame
            The DataFrame containing the CenSat BED data.
        sat_type : str
            The satellite type to filter the CenSat data.

        Returns:
        --------
        DataFrame
            A DataFrame containing the filtered CenSat data.

        Raises:
        -------
        ValueError
            If the filtering results in an empty DataFrame.
        '''
        if len(self.cenSat.columns) > 3:
            self.filtered_cenSat = cenSat[cenSat[3].str.contains(sat_type)]
            if self.filtered_cenSat.empty:
                raise ValueError("Filtering CenSat by the specified satellite type resulted in an empty DataFrame.")
        else:
            print('Regions File has no name column. Not filtering by sat_type.')
            self.filtered_cenSat = cenSat

        return self.filtered_cenSat


    def intersect_files(self, bed4Methyl, filtered_cenSat):
        '''
        Intersects the filtered bedMethyl and CenSat DataFrames using pybedtools.

        Parameters:
        -----------
        bedgraphMethyl : DataFrame
            The DataFrame containing the filtered bedMethyl data.
        subset_cenSat : DataFrame
            The DataFrame containing the filtered CenSat data.

        Returns:
        --------
        DataFrame
            A DataFrame containing the intersected regions between the filtered 
            bedMethyl and CenSat files.

        Raises:
        -------
        ValueError
            If the intersection results in an empty DataFrame.
        '''
        bedMethyl_bedtool = pybedtools.BedTool.from_dataframe(bed4Methyl)
        cenSat_bedtool = pybedtools.BedTool.from_dataframe(filtered_cenSat)
        intersected_bedtool = bedMethyl_bedtool.intersect(cenSat_bedtool, wa=True, u=True)
        self.intersected_bed4Methyl = intersected_bedtool.to_dataframe()
        
        if self.intersected_bed4Methyl.empty:
            raise ValueError("The intersection resulted in an empty DataFrame.")
        
        return self.intersected_bed4Methyl
    

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
    
    hmmCDRparser.read_cenSat(path=hmmCDRparser.cenSat_path)
    hmmCDRparser.read_bedMethyl(path=hmmCDRparser.bedMethyl_path)
    hmmCDRparser.filter_cenSat(cenSat=hmmCDRparser.cenSat, 
                               sat_type=hmmCDRparser.sat_type)
    hmmCDRparser.filter_bedMethyl(bedMethyl=hmmCDRparser.bedMethyl,
                                  mod_code=hmmCDRparser.mod_code)
    hmmCDRparser.intersect_files(bed4Methyl=hmmCDRparser.bed4Methyl,
                                 filtered_cenSat=hmmCDRparser.filtered_cenSat)
    
    # Save the filtered DataFrames to TSV files
    hmmCDRparser.intersected_bed4Methyl.to_csv(f'{args.output_prefix}_filtered_bedMethyl.bedgraph', sep='\t', index=False, header=False)
    hmmCDRparser.filtered_cenSat.to_csv(f'{args.output_prefix}_filtered_cenSat.bed', sep='\t', index=False, header=False)
    print(f"Filtered bedgraphMethyl saved to: {args.output_prefix}_filtered_bedgraphMethyl.bedgraph")
    print(f"Filtered CenSat saved to: {args.output_prefix}_filtered_cenSat.bed")