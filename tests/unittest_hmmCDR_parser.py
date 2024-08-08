import unittest
from unittest.mock import patch
import concurrent.futures
import pandas as pd
import tempfile
import os

from hmmCDRparse import hmmCDR_parser

class TestHMMCDRParser(unittest.TestCase):
    
    def setUp(self):
        # Create temporary files for the tests
        self.temp_bedMethyl = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.bed')
        self.temp_cenSat = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.bed')
        self.output_prefix = 'test_output'
    
    def tearDown(self):
        # Cleanup the temporary files after each test
        os.remove(self.temp_bedMethyl.name)
        os.remove(self.temp_cenSat.name)

    def test_bedMethyl_file_empty(self):
        """
        Test handling of an empty bedMethyl file.
        """
        # Writing an empty DataFrame to the temp file
        pd.DataFrame().to_csv(self.temp_bedMethyl.name, sep='\t', index=False, header=False)
        
        with self.assertRaises(ValueError):
            parser = hmmCDR_parser(bedMethyl_path=self.temp_bedMethyl.name,
                                    cenSat_path=self.temp_cenSat.name,
                                    mod_code='m',
                                    sat_type='H1L',
                                    bedgraph=False)
            parser.read_bedMethyl(path=parser.bedMethyl_path)

    def test_cenSat_file_empty(self):
        """
        Test handling of an empty cenSat file.
        """
        # Writing an empty DataFrame to the temp file
        pd.DataFrame().to_csv(self.temp_cenSat.name, sep='\t', index=False, header=False)
        
        with self.assertRaises(ValueError):
            parser = hmmCDR_parser(bedMethyl_path=self.temp_bedMethyl.name,
                                    cenSat_path=self.temp_cenSat.name,
                                    mod_code='m',
                                    sat_type='H1L',
                                    bedgraph=False)
            parser.read_cenSat(path=parser.cenSat_path)

    def test_filter_bedMethyl_valid(self):
        """
        Test filtering bedMethyl for a valid modification code.
        """
        # Writing valid bedMethyl data to the temp file
        bedMethyl_data = pd.DataFrame({
            0: ['chr1', 'chr1'],
            1: [100, 200],
            2: [150, 250],
            3: ['m', 'm'],
            4: ['m', 'm'],
            5: [0, 0], 6: [0, 0],
            7: [0, 0], 8: [0, 0], 9: [0, 0],
            10: [0.8, 0.9],
            11: [0, 0], 12: [0, 0], 13: [0, 0], 14: [0, 0], 
            15: [0, 0], 16: [0, 0], 17: [0, 0]
        })

        bedMethyl_data.to_csv(self.temp_bedMethyl.name, sep='\t', index=False, header=False)

        parser = hmmCDR_parser(bedMethyl_path=self.temp_bedMethyl.name,
                               cenSat_path='b',
                               mod_code='m',
                               sat_type='H1L',
                               bedgraph=False)
        
        bedMethyl = parser.read_bedMethyl(path=parser.bedMethyl_path)
        filtered_bedMethyl = parser.filter_bedMethyl(bedMethyl=bedMethyl, mod_code='m')
        self.assertFalse(filtered_bedMethyl.empty, "Filtered bedMethyl DataFrame should not be empty.")

    def test_filter_bedMethyl_valid(self):
        """
        Test reading and filtering bedMethyl for a valid modification code with the bedgraph flag
        """
        # Writing valid bedMethyl data to the temp file
        bedMethyl_data = pd.DataFrame({
            0: ['chr1', 'chr1'],
            1: [100, 200],
            2: [150, 250],
            10: [0.8, 0.9],
        })

        bedMethyl_data.to_csv(self.temp_bedMethyl.name, sep='\t', index=False, header=False)

        parser = hmmCDR_parser(bedMethyl_path=self.temp_bedMethyl.name,
                               cenSat_path='b',
                               mod_code='m',
                               sat_type='H1L',
                               bedgraph=True)
        
        bedMethyl = parser.read_bedMethyl(path=parser.bedMethyl_path)
        filtered_bedMethyl = parser.filter_bedMethyl(bedMethyl=bedMethyl, mod_code='m')
        self.assertFalse(filtered_bedMethyl.empty, "Filtered bedMethyl DataFrame should not be empty.")

    def test_filter_bedMethyl_invalid(self):
        """
        Test filtering bedMethyl with a modification code that doesn't exist in the file.
        """
        # Writing bedMethyl data without the modification code to the temp file
        bedMethyl_data = pd.DataFrame({
            0: ['chr1', 'chr1'],
            1: [100, 200],
            2: [150, 250],
            3: ['n', 'n'],
            10: [0.8, 0.9]
        })

        bedMethyl_data.to_csv(self.temp_bedMethyl.name, sep='\t', index=False, header=False)

        parser = hmmCDR_parser(bedMethyl_path=self.temp_bedMethyl.name,
                               cenSat_path='b',
                               mod_code='m',
                               sat_type='H1L',
                               bedgraph=False)
        
        bedMethyl = parser.read_bedMethyl(path=parser.bedMethyl_path)
        with self.assertRaises(ValueError):
            parser.filter_bedMethyl(bedMethyl=bedMethyl, mod_code='m')

    def test_filter_cenSat_valid(self):
        """
        Test filtering cenSat for a valid satellite type.
        """
        # Writing valid cenSat data to the temp file
        cenSat_data = pd.DataFrame({
            0: ['chr1', 'chr1'],
            1: [1000, 2000],
            2: [1500, 2500],
            3: ['H1L', 'H1L']
        })
        cenSat_data.to_csv(self.temp_cenSat.name, sep='\t', index=False, header=False)

        parser = hmmCDR_parser(bedMethyl_path='a',
                               cenSat_path=self.temp_cenSat.name,
                               mod_code='m',
                               sat_type='H1L',
                               bedgraph=False)
        
        cenSat = parser.read_cenSat(path=parser.cenSat_path)
        filtered_cenSat = parser.filter_cenSat(cenSat=cenSat, sat_type='H1L')
        self.assertFalse(filtered_cenSat.empty, "Filtered cenSat DataFrame should not be empty.")

    def test_filter_cenSat_invalid(self):
        """
        Test filtering cenSat with a satellite type that doesn't exist in the file.
        """
        # Writing cenSat data without the target satellite type to the temp file
        cenSat_data = pd.DataFrame({
            0: ['chr1', 'chr1'],
            1: [1000, 2000],
            2: [1500, 2500],
            3: ['XYZ', 'XYZ']
        })
        cenSat_data.to_csv(self.temp_cenSat.name, sep='\t', index=False, header=False)

        parser = hmmCDR_parser(bedMethyl_path=self.temp_bedMethyl.name,
                               cenSat_path=self.temp_cenSat.name,
                               mod_code='m',
                               sat_type='H1L',
                               bedgraph=False)
        
        cenSat = parser.read_cenSat(path=parser.cenSat_path)
        with self.assertRaises(ValueError):
            parser.filter_cenSat(cenSat=cenSat, sat_type='H1L')

    def test_intersect_files_valid(self):
        """
        Test intersecting valid bedgraphMethyl and cenSat subsets.
        """
        # Writing valid bedMethyl and cenSat data with overlapping regions to the temp files
        bedMethyl_data = pd.DataFrame({
            0: ['chr1', 'chr1'],
            1: [100, 101],
            2: [300, 301],
            3: ['m', 'm'],
            4: ['m', 'm'],
            5: [0, 0], 6: [0, 0],
            7: [0, 0], 8: [0, 0], 9: [0, 0],
            10: [0.8, 0.9],
            11: [0, 0], 12: [0, 0], 13: [0, 0], 14: [0, 0], 
            15: [0, 0], 16: [0, 0], 17: [0, 0]
        })
        cenSat_data = pd.DataFrame({
            0: ['chr1', 'chr1'],
            1: [150, 350],
            2: [250, 450],
            3: ['H1L', 'H1L']
        })
        bedMethyl_data.to_csv(self.temp_bedMethyl.name, sep='\t', index=False, header=False)
        cenSat_data.to_csv(self.temp_cenSat.name, sep='\t', index=False, header=False)

        parser = hmmCDR_parser(bedMethyl_path=self.temp_bedMethyl.name,
                               cenSat_path=self.temp_cenSat.name,
                               mod_code='m',
                               sat_type='H1L',
                               bedgraph=False)
        
        cenSat = parser.read_cenSat(path=parser.cenSat_path)
        bedMethyl = parser.read_bedMethyl(path=parser.bedMethyl_path)
        filtered_cenSat = parser.filter_cenSat(cenSat=cenSat, sat_type=parser.sat_type)
        filtered_bedMethyl = parser.filter_bedMethyl(bedMethyl=bedMethyl, mod_code=parser.mod_code)

        intersected_bed4Methyl = parser.intersect_files(bed4Methyl=filtered_bedMethyl, filtered_cenSat=filtered_cenSat)
        self.assertFalse(intersected_bed4Methyl.empty, "Intersected DataFrame should not be empty.")

    def test_intersect_files_no_overlap(self):
        """
        Test intersecting bedgraphMethyl and cenSat subsets with no overlap.
        """
        # Writing bedMethyl and cenSat data with non-overlapping regions to the temp files
        bedMethyl_data = pd.DataFrame({
            0: ['chr1', 'chr1'],
            1: [100, 500],
            2: [101, 501],
            3: ['m', 'm'],
            4: ['m', 'm'],
            5: [0, 0], 6: [0, 0],
            7: [0, 0], 8: [0, 0], 9: [0, 0],
            10: [0.8, 0.9],
            11: [0, 0], 12: [0, 0], 13: [0, 0], 14: [0, 0], 
            15: [0, 0], 16: [0, 0], 17: [0, 0]
        })
        cenSat_data = pd.DataFrame({
            0: ['chr1', 'chr1'],
            1: [150, 350],
            2: [250, 450],
            3: ['H1L', 'H1L']
        })
        bedMethyl_data.to_csv(self.temp_bedMethyl.name, sep='\t', index=False, header=False)
        cenSat_data.to_csv(self.temp_cenSat.name, sep='\t', index=False, header=False)

        parser = hmmCDR_parser(bedMethyl_path=self.temp_bedMethyl.name,
                               cenSat_path=self.temp_cenSat.name,
                               mod_code='m',
                               sat_type='H1L',
                               bedgraph=False)
        
        cenSat = parser.read_cenSat(path=parser.cenSat_path)
        bedMethyl = parser.read_bedMethyl(path=parser.bedMethyl_path)
        filtered_cenSat = parser.filter_cenSat(cenSat=cenSat, sat_type=parser.sat_type)
        filterd_bedMethyl = parser.filter_bedMethyl(bedMethyl=bedMethyl, mod_code=parser.mod_code)
        with self.assertRaises(ValueError):
            parser.intersect_files(bed4Methyl=filterd_bedMethyl, filtered_cenSat=filtered_cenSat)

    def test_parallel_processing(self):
        """
        Test that processing chromosomes in parallel produces the correct dictionary 
        structure and values.
        """
        # Create test data for bedMethyl and cenSat
        bedMethyl_data = pd.DataFrame({
            0: ['chr1', 'chr2', 'chr1', 'chr2'],
            1: [100, 200, 150, 250],
            2: [150, 250, 200, 300],
            3: ['m', 'm', 'm', 'm'],
            4: [0, 0, 0, 0],
            5: [0, 0, 0, 0],
            6: [0, 0, 0, 0],
            7: [0, 0, 0, 0],
            8: [0, 0, 0, 0],
            9: [0, 0, 0, 0],
            10: [0.8, 0.9, 0.85, 0.95],
            11: [0, 0, 0, 0],
            12: [0, 0, 0, 0],
            13: [0, 0, 0, 0],
            14: [0, 0, 0, 0],
            15: [0, 0, 0, 0],
            16: [0, 0, 0, 0],
            17: [0, 0, 0, 0]
        })

        cenSat_data = pd.DataFrame({
            0: ['chr1', 'chr2', 'chr1', 'chr2'],
            1: [140, 210, 190, 260],
            2: [160, 230, 210, 280],
            3: ['H1L', 'H1L', 'H1L', 'H1L']
        })

        bedMethyl_data.to_csv(self.temp_bedMethyl.name, sep='\t', index=False, header=False)
        cenSat_data.to_csv(self.temp_cenSat.name, sep='\t', index=False, header=False)

        # Initialize parser with test paths and arguments
        parser = hmmCDR_parser(bedMethyl_path=self.temp_bedMethyl.name,
                               cenSat_path=self.temp_cenSat.name,
                            mod_code='m',
                            sat_type='H1L',
                            bedgraph=False)
        
        bedMethyl = parser.read_bedMethyl(path=parser.bedMethyl_path)
        cenSat = parser.read_cenSat(path=parser.cenSat_path)
        
        # Call the process_all_chromosomes method
        bed4Methyl_chrom_dict, cenSat_chrom_dict = parser.process_all_chromosomes(bedMethyl=bedMethyl, cenSat=cenSat)

        # Verify that the dictionary contains the correct keys (chromosomes)
        self.assertIn('chr1', bed4Methyl_chrom_dict)
        self.assertIn('chr2', bed4Methyl_chrom_dict)
        self.assertIn('chr1', cenSat_chrom_dict)
        self.assertIn('chr2', cenSat_chrom_dict)

        # Verify the content of the intersected DataFrames for each chromosome
        self.assertFalse(bed4Methyl_chrom_dict['chr1'].empty)
        self.assertFalse(bed4Methyl_chrom_dict['chr2'].empty)
        self.assertFalse(cenSat_chrom_dict['chr1'].empty)
        self.assertFalse(cenSat_chrom_dict['chr2'].empty)


if __name__ == '__main__':
    unittest.main()