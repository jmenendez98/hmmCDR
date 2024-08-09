# hmmCDR


## Installation Methods: (need to add)
1. pyPI
2. conda
3. docker


This software is designed to find Centromere Dip Regions (CDRs), subCDRs, and their boundaries within the centromeric active alpha satellite (alpha-sat) array. CDRs are a uniquely hypo-methylated region within the typically hyper-methylated alpha-sat array. CDRs are tightly associated with the histone mark Centromere Protein A (CENP-A). This makes establishing accurate boundaries to CDRs and subCDRs essential to studying their relationship with CENPA. This method combines previous methods of identifying CDRs, through a sliding-window approach, with a Hidden Markov Model(HMM) that uses these sliding window estimates as a prior. The advantage to this two-fold approach is seen at the edges of the CDRs. A sliding window algorithm has a hard time drawing precise boundaries and identifying transitions in/out of the CDRs, whereas the HMM greatly improves identification of these regions. 

[include a photo of HMM improvement over sliding window]

[Include the photo I shared with Karen summarizing the HMM]

This python package takes in a bed file of 5mC methylation in aggregate, preferably from [modkit](https://github.com/nanoporetech/modkit), and an [Centromere-Satellite Annotation](https://github.com/kmiga/alphaAnnotation)(CenSat) file. The aggregate methylation file is used to determine where the 5mC depleted regions are, and the CenSat Annotation is used to subset the methylation files to only the alpha-sat array. This improves both the speed and accuracy of the CDR identification, as outside this region the trend of hypermethylation is not as strong. This package also processes each chromosome separately and in parallel to further improve speed.


## Inputs:
### 1. Modkit Pileup bedMethyl File:   

| column | name                  | description                                                                    | type  |
|--------|-----------------------|--------------------------------------------------------------------------------|-------|
| 1      | chrom                 | name of chromosome/contig                                                      | str   |
| 2      | start position        | 0-based start position                                                         | int   |
| 3      | end position          | 0-based exclusive end position                                                 | int   |
| 4      | modified base code    | single letter code for modified base                                           | str   |
| 5      | score                 | Equal to N<sub>valid_cov</sub>.                                                | int   |
| 6      | strand                | '+' for positive strand '-' for negative strand, '.' when strands are combined | str   |
| 7      | start position        | included for compatibility                                                     | int   |
| 8      | end position          | included for compatibility                                                     | int   |
| 9      | color                 | included for compatibility, always 255,0,0                                     | str   |
| 10     | N<sub>valid_cov</sub> | Refer to [modkit github](https://github.com/nanoporetech/modkit)               | int   |
| 11     | fraction modified     | N<sub>mod</sub> / N<sub>valid_cov</sub>                                        | float |
| 12     | N<sub>mod</sub>       | Refer to [modkit github](https://github.com/nanoporetech/modkit)               | int   |
| 13     | N<sub>canonical</sub> | Refer to [modkit github](https://github.com/nanoporetech/modkit)               | int   |
| 14     | N<sub>other_mod</sub> | Refer to [modkit github](https://github.com/nanoporetech/modkit)               | int   |
| 15     | N<sub>delete</sub>    | Refer to [modkit github](https://github.com/nanoporetech/modkit)               | int   |
| 16     | N<sub>fail</sub>      | Refer to [modkit github](https://github.com/nanoporetech/modkit)               | int   |
| 17     | N<sub>diff</sub>      | Refer to [modkit github](https://github.com/nanoporetech/modkit)               | int   |
| 18     | N<sub>nocall</sub>    | Refer to [modkit github](https://github.com/nanoporetech/modkit)               | int   |

### 2.  CenSat Annotation bed

| column | name                  | description                                                                    | type  |
|--------|-----------------------|--------------------------------------------------------------------------------|-------|
| 1      | chrom                 | name of chromosome/contig                                                      | str   |
| 2      | start position        | 0-based start position                                                         | int   |
| 3      | end position          | 0-based exclusive end position                                                 | int   |
| 4      | satellite type/name   | type of satellite and for some specific name in parentheses                    | str   |
| 5      | score                 | Not sure what if it is used for anytime.                                       | int   |
| 6      | strand                | '+' for positive strand '-' for negative strand, '.' if uncertain              | str   |
| 7      | start position        | included for compatibility                                                     | int   |
| 8      | end position          | included for compatibility                                                     | int   |
| 9      | color                 | color of the annotation in browser                                             | str   |


### Help Documentation
```
usage: hmmCDR [-h] [--mod_code MOD_CODE] [--sat_type SAT_TYPE]
              [--min_valid_cov MIN_VALID_COV] [--bedgraph]
              [--window_size WINDOW_SIZE]
              [--priorCDR_percent PRIORCDR_PERCENT]
              [--priorTransition_percent PRIORTRANSITION_PERCENT]
              [--minCDR_size MINCDR_SIZE] [--enrichment]
              [--cdr_priors CDR_PRIORS | --emission_matrix EMISSION_MATRIX | --transition_matrix TRANSITION_MATRIX]
              [--use_percentiles] [--n_iter N_ITER] [-w W] [-x X] [-y Y]
              [-z Z] [--save_intermediates] [--output_label OUTPUT_LABEL]
              bedMethyl_path cenSat_path output_path

Process input files with optional parameters.

positional arguments:
  bedMethyl_path        Path to the bedMethyl file
  cenSat_path           Path to the CenSat BED file
  output_path           Output Path for the output files

optional arguments:
  -h, --help            show this help message and exit
  --mod_code MOD_CODE   Modification code to filter bedMethyl file. (default: "m")
  --sat_type SAT_TYPE   Satellite type/name to filter CenSat bed file. (default: "H1L")
  --min_valid_cov MIN_VALID_COV Minimum Valid Coverage to consider a methylation site. (default: 10)
  --bedgraph            Flag indicating if the input is a bedgraph. (default: False)
  --window_size WINDOW_SIZE Window size to calculate prior regions. (default: 1020)
  --priorCDR_percent PRIORCDR_PERCENT   Percentile for finding priorCDR regions. (default: 5)
  --priorTransition_percent PRIORTRANSITION_PERCENT Percentile for finding priorTransition regions. (default: 10)
  --minCDR_size MINCDR_SIZE Minimum size for CDR regions. (default: 3000)
  --enrichment          Enrichment flag. Pass in if you are looking for methylation enriched regions. (default: False)
  --cdr_priors CDR_PRIORS   Path to the priorCDR bedfile
  --emission_matrix EMISSION_MATRIX Path to the emission matrix TSV file
  --transition_matrix TRANSITION_MATRIX Path to the transition matrix TSV file
  --use_percentiles     Use values for flags w,x,y,z as percentile cutoffs for each category. (default: False)
  --n_iter N_ITER       Maximum number of iteration allowed for the HMM. (default: 1)
  -w W                  Theshold for methylation to be classified as very low (default: 0)
  -x X                  Theshold for methylation to be classified as low (default: 25)
  -y Y                  Theshold for methylation to be classified as medium (default: 50)
  -z Z                  Theshold for methylation to be classified as high (default: 75)
  --save_intermediates  Set to true if you would like to save intermediates(filtered beds+window means). (default: False)
  --output_label OUTPUT_LABEL   Label to use for name column of hmmCDR BED file. Needs to match priorCDR label. (default: "CDR")
```

```
usage: hmmCDRprior [-h] [-m MOD_CODE] [-s SAT_TYPE] [--bedgraph]
                   [-w WINDOW_SIZE] [--priorCDR_percent PRIORCDR_PERCENT]
                   [--priorTransition_percent PRIORTRANSITION_PERCENT]
                   [--minCDR_size MINCDR_SIZE] [--enrichment]
                   [--save_intermediates] [--output_label OUTPUT_LABEL]
                   bedMethyl_path cenSat_path output_path

Process bedMethyl and CenSat BED file to produce hmmCDR priors

positional arguments:
  bedMethyl_path        Path to the bedMethyl file
  cenSat_path           Path to the CenSat BED file
  output_path           Path to the output priorCDRs BED file

optional arguments:
  -h, --help            show this help message and exit
  -m MOD_CODE, --mod_code MOD_CODE  Modification code to filter bedMethyl file (default: "m")
  -s SAT_TYPE, --sat_type SAT_TYPE  Satellite type/name to filter CenSat bed file. (default: "H1L")
  --bedgraph            Flag indicating if the input is a bedgraph. (default: False)
  -w WINDOW_SIZE, --window_size WINDOW_SIZE Window size to calculate prior regions. (default: 1020)
  --priorCDR_percent PRIORCDR_PERCENT   Percentile for finding priorCDR regions. (default: 5)
  --priorTransition_percent PRIORTRANSITION_PERCENT Percentile for finding priorTransition regions. (default: 10)
  --minCDR_size MINCDR_SIZE Minimum size for CDR regions. (default: 3000)
  --enrichment          Enrichment flag. Pass in if you are looking for methylation enriched regions. (default: False)
  --save_intermediates  Set to true if you would like to save intermediates(filtered beds+window means). (default: False)
  --output_label OUTPUT_LABEL Label to use for name column of priorCDR BED file. (default: "CDR")
```

```
usage: hmmCDRparse [-h] [--bedgraph] [--min_valid_cov MIN_VALID_COV]
                   [-m MOD_CODE] [-s SAT_TYPE]
                   bedMethyl_path cenSat_path output_prefix

Process bedMethyl and CenSat BED file to produce hmmCDR priors

positional arguments:
  bedMethyl_path        Path to the bedMethyl file
  cenSat_path           Path to the CenSat BED file
  output_prefix         Path to the output priorCDRs BED file

optional arguments:
  -h, --help            show this help message and exit
  --bedgraph            Flag indicating if the input is a bedgraph. (default: False)
  --min_valid_cov MIN_VALID_COV Minimum Valid Coverage to consider a methylation site. (default: 10)
  -m MOD_CODE, --mod_code MOD_CODE  Modification code to filter bedMethyl file (default: "m")
  -s SAT_TYPE, --sat_type SAT_TYPE  Satellite type/name to filter CenSat bed file. (default: "H1L")
```