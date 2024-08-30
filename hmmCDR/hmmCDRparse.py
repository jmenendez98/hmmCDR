import pandas as pd
import pybedtools
import os


class hmmCDRparse:
    def __init__(self, min_valid_cov, mod_code, sat_type, bedgraph=False):
        self.mod_code = mod_code
        self.sat_type = sat_type
        self.bedgraph = bedgraph
        self.min_valid_cov = min_valid_cov

    def read_bedMethyl(self, path):
        if not os.path.exists(path) or os.stat(path).st_size == 0:
            raise ValueError(f"bedMethyl file {path} does not exist or is empty.")
        
        bedMethyl = pd.read_csv(path, sep='\t', header=None, index_col=None)
        expected_columns = 4 if self.bedgraph else 18
        if len(bedMethyl.columns) != expected_columns:
            raise ValueError(f"Valid {'bedgraph' if self.bedgraph else 'bedMethyl'} should have {expected_columns} columns.")

        if self.bedgraph:
            return bedMethyl

        bed4Methyl = bedMethyl[(bedMethyl[3] == self.mod_code) & 
                                (bedMethyl[4] >= self.min_valid_cov)].iloc[:, [0, 1, 2, 10]]

        if bed4Methyl.empty:
            raise ValueError("Filtering bedMethyl by the specified modification code resulted in an empty DataFrame.")
        
        bed4Methyl.columns = ['chrom', 'start', 'end', 'name']

        return bed4Methyl
    
    def read_cenSat(self, path):
        if not os.path.exists(path) or os.stat(path).st_size == 0:
            raise ValueError(f"cenSat file {path} does not exist or is empty.")
        
        cenSat = pd.read_csv(path, sep='\t', header=None, index_col=None)

        if len(cenSat.columns) > 3:
            cenSat = cenSat[cenSat[3].str.contains('|'.join(self.sat_type))]
            if cenSat.empty:
                raise ValueError("Filtering CenSat by the specified satellite type resulted in an empty DataFrame.")
        else:
            print('Regions File has no name column. Not filtering by sat_type.')

        cenSat = cenSat.iloc[:, :4]
        cenSat.columns = ['chrom', 'start', 'end', 'name']

        return cenSat

    def intersect_files(self, bed4Methyl, cenSat):
        bedMethyl_bedtool = pybedtools.BedTool.from_dataframe(bed4Methyl)
        cenSat_bedtool = pybedtools.BedTool.from_dataframe(cenSat)
        intersected_bed4Methyl = bedMethyl_bedtool.intersect(cenSat_bedtool, wa=True, u=True).to_dataframe()
        
        if intersected_bed4Methyl.empty:
            raise ValueError("The intersection resulted in an empty DataFrame.")
        
        return intersected_bed4Methyl
    
    def sep_by_chrom(self, intersected_bed4Methyl, cenSat):
        chromosomes = cenSat['chrom'].unique()
        return (
            {chrom: intersected_bed4Methyl[intersected_bed4Methyl['chrom'] == chrom].reset_index(drop=True) for chrom in chromosomes},
            {chrom: cenSat[cenSat['chrom'] == chrom].reset_index(drop=True) for chrom in chromosomes}
        )

    def process_files(self, bedMethyl_path, cenSat_path):
        bed4Methyl = self.read_bedMethyl(bedMethyl_path)
        cenSat = self.read_cenSat(cenSat_path)
        intersected_bed4Methyl = self.intersect_files(bed4Methyl, cenSat)

        return self.sep_by_chrom(intersected_bed4Methyl, cenSat)