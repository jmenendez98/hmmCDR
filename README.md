# CDR_workflow
This is a WDL workflow to identify Centromere Dip Regions(CDRs)/Centromere Cores within the Active Alpha Sat Array. It does this by leveraging ONT's [modkit](https://github.com/nanoporetech/modkit) software to aggregate 5mC methylation encoded in MM/ML tags. This workflow and software was designed for use with HG002 Centomere Projects and the Human Pangenome Research Consortium's assemblies, to provide a highly accurate prediction for the CDR location. 

#### WDL Workflow Inputs:
    Reference Fasta: Reference used for the BAM files alignment.    
    Censat Bed: CenSat Track.    
    5mC Bam: BAM file containing 5mC tags in MM/ML form.     
    Sample ID: Sample ID, used for output/intermediate file naming convention.    

#### WDL Steps:
    1. Extract H1L
    2. Index BAM File
    3. Modkit Pileup
    4. Scatter on Contig Names
        5. Extract Individual Contigs from methylBed + bedtools intersect with H1L
        6. Strict CDR Detection
        7. HMM CDR Detection
        8. Create Validation Plots
    9. Gather Contigs and create outputs

#### Outputs:
    HMM CDR Bed: Output of workflow, HMM CDRs in bed9 file.
    methylBed: Whole genome bedMethyl File.
    
    IGV tar.gz: Contains H1L intersected bedMethyls and Strict CDRs that HMM was trained on.
    Summary tag.gz: Contains png files for each contig of validation track pngs, heatmaps, and histograms.

## Scripts:
#### strict_scoring.sh    
Thank you [Mira Mastoras](https://github.com/miramastoras) for starting this automated [CDR Detection](https://github.com/miramastoras/CDR_detect/tree/main) work! This script processes a bed4 of the fraction_modified generated from a modkit pileup and generates a CDR prediction file based on scoring of windows. Scoring happens at each percentile in the range provided. Windows with the highest scores are labeled as CDRs and transitions. (Note: Requires bedtools)
```
usage: strict_score.sh [-h]
    -i, --input: Input bed4 containing fraction_modified information from modkit pileup. [required]
    -r, --hg002_merged_H1L: Bed Track of regions you want to find CDRs within. CDRs are identified by regions of aggregate depleted 5mC methylation [required]
    -o, --output_prefix: A prefix for the output files. The output file will be named as <output_prefix>.strictScores.bed. [Required]
    --low_percentage: The starting percentile for scoring windows. [default 1]
    --high_percentage: Optional. The ending percentile for scoring windows. [default 20]
    --high-confidence: Cutoff for high confidence CDR score. Scores equal or above are high confidence. [default 17]
    --med-confidence: Cutoff for medium confidence CDR score. Scores above are medium confidence. [default 10]
    --low-confidence: Cutoff for low confidence CDR score. Scores above are low confidence. [default 5]
```

#### HMMCDRDetection.py
Thank you [Justin Chan](https://github.com/Justinmchan408) for writing the framework for this [HMM](https://github.com/Justinmchan408/HMMCDRDetection)! This version contains a few bug fixes, performance changes, and functionality improvements. One of the major ones is the ability to identify CDR transtions which are regions flanking CDRs where methylation gradually increases or decreases. (Note: Requires python packages numpy and pandas)
```
usage: HMMCDRReferenceDetection.py [-h]  
    -p bed4 file containing modified CPG site probabilities [required]
    -s bed file containing estimate CDR and Transition Regions [required]
    -o output bed prefix [required]
    -l Learning Rate for the Viterbi Learning [default 0.0000001]
    --steps Maximum steps for Viterbi Learning [default 100]
```

#### cdr_histogram.py
Creates a histogram of the `modkit pileup` values within the H1L array. Colors values based on whether they fall within annotations from the HMM CDR bed9 file. 
```
usage: cdr_histogram.py [-h]  
    -i bed4 file containing modified CPG site probabilities [required]
    -r bed file containing HMM CDR and Transition Region Predictions [required]
    -o Output file [default 'histogram.png']
```

#### hmm_heatmaps.py
Generates two heatmaps to represent the emission and transition matrices from HMM CDR predictions. 
```
usage: hmm_heatmaps.py [-h]  
    -e Emmission matrix in .csv format from HMMCDRDetection.py [required]
    -t Transition matrix in .csv format from HMMCDRDetection.py [required]
    -o Output prefix [default 'hmm_heatmap']
```

