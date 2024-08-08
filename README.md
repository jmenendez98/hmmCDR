# hmmCDR

Input Requirements:
* Modkit Pileup bedMethyl
* CenSat Annotation bed


# WIP:

Testing Commands:
```
# testing the parser
python3 hmmCDR_parser.py --mod_code m --sat_type H1L chr10_MAT_HG002_ONT.5mCpileup.bed chr10_MAT_hg002v1.0.1.cenSatv2.0.bed chr10_MAT_hmmCDR_parser_test
# testing the priors
python3 hmmCDR_priors.py --mod_code m --sat_type H1L chr10_MAT_HG002_ONT.5mCpileup.bed chr10_MAT_hg002v1.0.1.cenSatv2.0.bed chr10_MAT_hmmCDR_priors_test.bed
# w/ bedgraph input
python3 hmmCDR_priors.py --mod_code m --bedgraph --sat_type H1L chr10_hmmCDR_parser_test_filtered_bedMethyl.bedgraph chr10_MAT_hg002v1.0.1.cenSatv2.0.bed chr10_hmmCDR_priors_test_bedg.bed
# testing the actual HMM
python3 hmmCDR.py --mod_code m --sat_type H1L chr10_MAT_HG002_ONT.5mCpileup.bed chr10_MAT_hg002v1.0.1.cenSatv2.0.bed chr10_MAT_hmmCDR_full_test.bed
```

# What to do next:
1. DocStrings!!!
2. Make sure all the flags are working as intended
    - Can I run with a matrix input?
    - Can I run with a priorCDR input?
