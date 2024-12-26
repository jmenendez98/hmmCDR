import os

import pytest

from hmmCDR.bed_parser import bed_parser


class TestParser:
    @pytest.fixture
    def test_data_dir(self):
        return os.path.join("tests", "data")

    @pytest.fixture
    def parser(self):
        return bed_parser(
            mod_code="m",
            bedgraph=False,
            min_valid_cov=10,
            sat_type=["active_hor"],
            pre_subset_censat=False,
        )

    def test_fake_bedfile(self, test_data_dir, parser):
        """Test handling of nonexistent file"""
        nonexistent_file = os.path.join(test_data_dir, "nonexistent.bed")
        with pytest.raises(FileNotFoundError):
            parser.read_and_filter_bedfile(nonexistent_file, file_type="bedmethyl")
        with pytest.raises(FileNotFoundError):
            parser.read_and_filter_bedfile(nonexistent_file, file_type="censat")

    def test_empty_bedfile(self, test_data_dir, parser):
        """Test handling of empty file"""
        empty_file = os.path.join(test_data_dir, "empty.bed")

        result1 = parser.read_and_filter_bedfile(str(empty_file), file_type="bedmethyl")
        result2 = parser.read_and_filter_bedfile(str(empty_file), file_type="censat")

        assert len(result1) == 0
        assert len(result2) == 0

    def test_censat_bedfile(self, test_data_dir, parser):
        """Test basic censat reading functionality"""
        sample_censat_bed = os.path.join(test_data_dir, "censat_test.bed")
        results = parser.read_and_filter_bedfile(sample_censat_bed, file_type="censat")

        print(results)

        assert len(results) == 1
        assert results[0][0] == "chrX_MATERNAL"
        assert int(results[0][1]) == 57866525
        assert int(results[0][2]) == 60979767
        assert results[0][3] == "active_hor(S3CXH1L)"

    def test_bedmethyl_bedfile(self, test_data_dir, parser):
        """Test basic bedmethyl reading functionality"""
        sample_bedmethyl_bed = os.path.join(test_data_dir, "bedmethyl_test.bed")
        results = parser.read_and_filter_bedfile(
            sample_bedmethyl_bed, file_type="bedmethyl"
        )

        assert len(results) == 62064

    def test_chrom_dict_len(self, test_data_dir, parser):
        """Test basic bedmethyl reading functionality"""
        sample_bedmethyl_bed = os.path.join(test_data_dir, "bedmethyl_test.bed")
        sample_censat_bed = os.path.join(test_data_dir, "censat_test.bed")

        methylation_chrom_dict, regions_chrom_dict = parser.process_files(
            bedmethyl_path=sample_bedmethyl_bed,
            censat_path=sample_censat_bed,
        )

        # lengths should be 1 because test file is only one chromosome
        assert len(methylation_chrom_dict) == 1
        assert len(regions_chrom_dict) == 1
