# hmmCDR

This software is designed to find Centromere Dip Regions (CDRs), subCDRs, and their boundaries within the centromeric active alpha satellite (alpha-sat) array. CDRs are a uniquely hypo-methylated region within the typically hyper-methylated alpha-sat array. CDRs are tightly associated with the histone mark Centromere Protein A (CENP-A). This makes establishing accurate boundaries to CDRs and subCDRs essential to studying their relationship with CENPA. This method combines previous methods of identifying CDRs, through a sliding-window approach, with a Hidden Markov Model(HMM) that uses these sliding window estimates as a prior. The advantage to this two-fold approach is seen at the edges of the CDRs. A sliding window algorithm has a hard time drawing precise boundaries and identifying transitions in/out of the CDRs, whereas the HMM greatly improves identification of these regions. 

[include a photo of HMM improvement over sliding window]

This python package takes in a bed file of 5mC methylation in aggregate, preferably from [modkit](https://github.com/nanoporetech/modkit), and an [Centromere-Satellite Annotation](https://github.com/kmiga/alphaAnnotation)(CenSat) file. The aggregate methylation file is used to determine where the 5mC depleted regions are, and the CenSat Annotation is used to subset the methylation files to only the alpha-sat array. This improves both the speed and accuracy of the CDR identification, as outside this region the trend of hypermethylation is not as strong. This package also processes each chromosome separately and in parallel to further improve speed.

## Installation Methods:

## Prefered Inputs:
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
| 5      | score                 | Equal to N<sub>valid_cov</sub>.                                                | int   |
| 6      | strand                | '+' for positive strand '-' for negative strand, '.' if uncertain              | str   |
| 7      | start position        | included for compatibility                                                     | int   |
| 8      | end position          | included for compatibility                                                     | int   |
| 9      | color                 | color of the annotation in browser                                             | str   |

