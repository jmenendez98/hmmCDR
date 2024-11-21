# `hmmCDR`

[![Conda](https://img.shields.io/conda/vn/jmmenend/hmmcdr?label=conda&color=green)](https://anaconda.org/jmmenend/hmmcdr)
[![DockerHub](https://img.shields.io/docker/v/jmmenend/hmmcdr?label=DockerHub&color=blue)](https://hub.docker.com/r/jmmenend/hmmcdr)
[![pypi](https://img.shields.io/pypi/v/hmmCDR)](https://pypi.org/project/hmmCDR/0.1.4/)

`hmmCDR` is a set of python scripts to automate the prediction of hypomethylated Centromere Dip Regions (CDRs) within active alpha satellite HOR arrays. Utilizes an Hidden Markov Model to allow finer resolution of subCDR boundaries.

**Inputs:**      
1. bedmethyl, preferably generated with [`modkit`](https://github.com/nanoporetech/modkit) `pileup` (Or bedgraph of 5mC methylation with the `--bedgraph` flag)
2. Centromeric Satellite Annotations, designed for use with annotations generated using https://github.com/kmiga/alphaAnnotation.
3. Output file name

## Installation: 

`hmmCDR` can be installed through `conda`. With `bioconda` and `conda-forge` channels enabled.
```bash
conda install jmmenend::hmmcdr -c bioconda -c conda-forge
```

`hmmCDR` can be run with `docker`.
```bash
docker run -v .:/data jmmenend/hmmcdr:0.2.3 # INPUTS/FLAGS # 
```

`hmmCDR` can be install through `pypi`. **NOTE**: This requires a separate installation of `bedtools` in the environment.
```bash
pip install hmmCDR
```

## Description:

`<img src="imgs/CDR_HMM_diagram.png" alt="HMM Diagram" width="600">`
* **Replace with figure 1 from HMM CDR Paper...**              

This software is designed to find Centromere Dip Regions (CDRs), subCDRs, and their boundaries within the centromeric active alpha satellite (alpha-sat) array. CDRs are a uniquely hypo-methylated region within the typically hyper-methylated alpha-sat array. CDRs are tightly associated with the histone mark Centromere Protein A (CENP-A). This makes establishing accurate boundaries to CDRs and subCDRs essential to studying their relationship with CENPA. This method combines previous methods of identifying CDRs, through a sliding-window approach, with a Hidden Markov Model (HMM) that uses these sliding window estimates as a prior. The advantage to this two-fold approach is seen at the edges of the CDRs. A sliding window algorithm has a hard time drawing precise boundaries and identifying transitions in/out of the CDRs, whereas the HMM greatly improves identification of these regions. 

This python package takes in a bed file of 5mC methylation in aggregate, preferably from [modkit](https://github.com/nanoporetech/modkit), and an [Centromere-Satellite Annotation](https://github.com/kmiga/alphaAnnotation) (CenSat) file. The aggregate methylation file is used to determine where the 5mC depleted regions are, and the CenSat Annotation is used to subset the methylation files to only the alpha-sat array. This improves both the speed and accuracy of the CDR identification, as outside this region the trend of hypermethylation is not as strong. This package also processes each chromosome separately and in parallel to further improve speed.

## Input Details:
### 1. Modkit bedMethyl (Refer to https://github.com/nanoporetech/modkit)    

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

### 2.  CenSat Annotation (Refer to https://github.com/kmiga/alphaAnnotation)         

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


### Usage:
```
usage: hmmCDR [-h] [-m MOD_CODE] [--bedgraph] [--min_valid_cov MIN_VALID_COV] [-s SAT_TYPE] [--pre_subset_censat] [--window_size WINDOW_SIZE] [--step_size STEP_SIZE] [--prior_threshold PRIOR_THRESHOLD] [--prior_use_percentile] [--min_prior_size MIN_PRIOR_SIZE] [--enrichment]
              [--percentile_emissions] [-w W] [-x X] [-y Y] [-z Z] [--e_matrix E_MATRIX] [--t_matrix T_MATRIX] [--n_iter N_ITER] [--tol TOL] [--hmm_merge_distance HMM_MERGE_DISTANCE] [--min_cdr_size MIN_CDR_SIZE] [--min_cdr_score MIN_CDR_SCORE]
              [--min_low_conf_size MIN_LOW_CONF_SIZE] [--min_low_conf_score MIN_LOW_CONF_SCORE] [--main_color MAIN_COLOR] [--low_conf_color LOW_CONF_COLOR] [--min_subCDRs MIN_SUBCDRS] [--large_merge_distance LARGE_MERGE_DISTANCE] [--output_all] [--output_label OUTPUT_LABEL]
              bedmethyl censat output

Process input files with optional parameters.

positional arguments:
  bedmethyl             Path to the bedMethyl file
  censat                Path to the CenSat BED file
  output                Output Path for the output files

options:
  -h, --help            show this help message and exit
  -m MOD_CODE, --mod_code MOD_CODE
                        Modification code to filter bedMethyl file (default: "m")
  --bedgraph            Flag indicating if the input is a bedgraph. (default: False)
  --min_valid_cov MIN_VALID_COV
                        Minimum valid coverage to consider a methylation site (read from full modkit pileup files). (default: 10)
  -s SAT_TYPE, --sat_type SAT_TYPE
                        Comma-separated list of satellite types/names to filter CenSat bed file. (default: "H1L")
  --pre_subset_censat   Set flag if your annotations bed file is already subset to only the region you desire. (default: False)
  --window_size WINDOW_SIZE
                        Window size to calculate prior regions. (default: 1190)
  --step_size STEP_SIZE
                        Step size when calculation windows for priors. (default: 1190)
  --prior_threshold PRIOR_THRESHOLD
                        Threshold for determining if a window is a CDR. Uses this percentile if --prior_use_percentile is passed (default: 30.0)
  --prior_use_percentile
                        Whether or not to use percentile when calculating windowing priors. (default: False)
  --min_prior_size MIN_PRIOR_SIZE
                        Minimum size for CDR regions. (default: 8330)
  --enrichment          Enrichment flag. Pass in if you are looking for methylation enriched regions. (default: False)
  --percentile_emissions
                        Use values for flags w,x,y,z as raw threshold cutoffs for each emission category. (default: False)
  -w W                  Threshold of non-zero methylation percentile to be classified as None (default: 0.0)
  -x X                  Threshold of non-zero methylation percentile to be classified as low (default: 33.3)
  -y Y                  Threshold of non-zero methylation percentile to be classified as medium (default: 66.6)
  -z Z                  Threshold of non-zero methylation percentile to be classified as high (default: 100.0)
  --e_matrix E_MATRIX   Custom Emission Matrix (Ex: [[0.002,0.10,0.28,0.60],[0.05,0.85,0.08,0.02]])
  --t_matrix T_MATRIX   Custom Transition Matrix (Ex: [[0.9999,0.003],[0.0001,0.997]])
  --n_iter N_ITER       Maximum number of iteration allowed for the HMM. (default: 1)
  --tol TOL             Cutoff for model convergence in hmmlearn. (default: 10)
  --hmm_merge_distance HMM_MERGE_DISTANCE
                        Distance to merge adjacently labelled subCDR regions. (default: 1190)
  --min_cdr_size MIN_CDR_SIZE
                        Minimum size of region identified. (default: 1190)
  --min_cdr_score MIN_CDR_SCORE
                        The minimum HMM score [0-100] required to call a CDR. (default: 95)
  --min_low_conf_size MIN_LOW_CONF_SIZE
                        Minimum size of region identified. (default: 0)
  --min_low_conf_score MIN_LOW_CONF_SCORE
                        The minimum HMM score [0-100] required to call a low confidence CDR. (default: 75)
  --main_color MAIN_COLOR
                        Color to dictate main regions. (default: 50,50,255)
  --low_conf_color LOW_CONF_COLOR
                        Color to dictate low confidence regions. (default: 100,150,200)
  --min_subCDRs MIN_SUBCDRS
                        Minimum number of subCDRs to report a CDR. (default: 3)
  --large_merge_distance LARGE_MERGE_DISTANCE
                        Distance to merge subCDRs into a larger CDR annotation. (default: 200000)
  --output_all          Set to true if you would like to save all intermediate filesf. (default: False)
  --output_label OUTPUT_LABEL
                        Label to use for name column of hmmCDR BED file. Needs to match priorCDR label. (default: "subCDR")
```

## License

This project is licensed under the MIT License.