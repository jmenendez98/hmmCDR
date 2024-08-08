# hmmCDR

Input Requirements:
* Modkit Pileup bedMethyl
* CenSat Annotation bed


# WIP:

Testing Commands:
```
# testing the parser
python3 hmmCDRparse.py --mod_code m --sat_type H1L ../tests/chrX_HG002_ONT.5mCpileup.bed ../tests/chrX_hg002v1.0.1.cenSatv2.0.bed chrX_hmmCDR_parser_test

# testing the priors
python3 hmmCDRpriors.py --mod_code m --sat_type H1L ../tests/chrX_HG002_ONT.5mCpileup.bed ../tests/chrX_hg002v1.0.1.cenSatv2.0.bed ../tests/chrX_hmmCDR_priors_test.bed
# w/ bedgraph input
python3 hmmCDR_priors.py --mod_code m --bedgraph --sat_type H1L chr10_hmmCDR_parser_test_filtered_bedMethyl.bedgraph chr10_MAT_hg002v1.0.1.cenSatv2.0.bed chr10_hmmCDR_priors_test_bedg.bed
# testing the actual HMM
python3 hmmCDR.py --mod_code m --sat_type H1L ../tests/chrX_HG002_ONT.5mCpileup.bed ../tests/chrX_hg002v1.0.1.cenSatv2.0.bed ../tests/chrX_hmmCDR_full_test.bed

### UNIT TESTING:
cd tests
python -m unittest -v unittest_hmmCDR_parser.py
python -m unittest -v unittest_hmmCDR_priors.py
python -m unittest -v unittest_hmmCDR.py
```

# What to do next:
1. DocStrings!!!
2. Make sure all the flags are working as intended
    -- 
    - Can I run with a priorCDR input?
3. Test with invalid inputs to throw decent errors:
    - multiple chromosomes in one bedfile
    - file not overlapping with region of interest
