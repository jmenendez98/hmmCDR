# CDR_workflow
This is a WDL workflow to identify Centromere Dip Regions(CDRs)/Centromere Cores within the Active Alpha Sat Array. It does this by leveraging ONT's [modkit](https://github.com/nanoporetech/modkit) software to aggregate 5mC methylation encoded in MM/ML tags. This workflow and software was designed for use with HG002 Centomere Projects and the Human Pangenome Research Consortium's assemblies, to provide a highly accurate prediction for the CDR location. 

#### WDL Workflow Inputs:
    Reference Fasta: Reference used for the BAM files alignment.    
    Censat Bed: CenSat Track.    
    5mC Bam: BAM file containing 5mC tags in MM/ML form.     
    Sample ID: Sample ID, string prefix used for output/intermediate file naming convention.    

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

## Running the WDL:
You can launch the WDL using files present in this repo. Using file in the [`wdl`](./wdl) folder. You can launch the CDR workflow with your own inputs, by editing the `inputs.json` for your data. Then, in an interactive session with at least 16 CPUs(For GI Private Cluster: `srun --job-name=interactive_medium --nodes=1 --cpus-per-task=16 --mem=256G --time=12:00:00 --partition=medium --pty bash`), run `bash cdr_detection.single_machine.wdl.sh`. (Assumes you have `toil` installed/configured properly refer to: [giwiki/toil](https://giwiki.gi.ucsc.edu/index.php?title=Phoenix_WDL_Tutorial))

## Scripts:

### CDR Detection:

#### strictCDRDetection.sh    
Thank you [Mira Mastoras](https://github.com/miramastoras) for starting this automated [CDR Detection](https://github.com/miramastoras/CDR_detect/tree/main) work! This script processes a bed4 of the fraction_modified generated from a modkit pileup and generates a CDR prediction file based on scoring of windows. Scoring happens at each percentile in the range provided. Windows with the highest scores are labeled as CDRs and transitions. (Note: Requires bedtools)
```
Usage: strict_cdr.sh [-h]
    -i, --input: Modkit pileup bedgraph of fraction modified [Required]
    -r, --hg002_merged_H1L: A BED file containing ONLY H1L(active-Alpha Satellite) regions [Required]
    -o, --output_prefix: A prefix for the output files. The output files will be named as <output_prefix>.strictCDR.bed and <output_prefix>.strictTransitions.bed [Required]
    -p, --percent: The percentage threshold for the CDRs. Default is 10 [Optional]
    -t, --transition_percent: The transition percentage threshold. Default is 20 [Optional]
    -w, --window_size: The window size. Default is 1020 [Optional]
    -m, --min_size: The minimum amound of windows flagged as a strict CDRs to keep a CDR. Minimum CDR Size is this * window_size. Default is 3 [Optional]
    -d, --merge_distance: The distance to merge nearby CDR entries(before filtering by size). Distance is given in number of windows. Default is 2 [Optional]
```

#### HMMCDRDetection.py
Thank you [Justin Chan](https://github.com/Justinmchan408) for writing the framework for this [HMM](https://github.com/Justinmchan408/HMMCDRDetection)! This version contains a few bug fixes, performance changes, and functionality improvements. One of the major ones is the ability to identify CDR transtions which are regions flanking CDRs where methylation gradually increases or decreases. (Note: Requires python packages numpy and pandas)
```
usage: HMMCDRReferenceDetection.py [-h]  
    -p: bed4 file containing modified CPG site probabilities [Required]
    -s: bed file containing estimate CDR and Transition Regions [Required]
    -o: output bed prefix [Required]
    -l: Learning Rate for the Viterbi Learning. Default is 0.0000001 [Optional]
    --steps: Maximum steps for Viterbi Learning. Default is 100 [Optional]
```

### HMM Validations:

#### cdr_histogram.py
Creates a histogram of the `modkit pileup` values within the H1L array. Colors values based on whether they fall within annotations from the HMM CDR bed9 file. 
```
usage: cdr_histogram.py [-h]  
    -i, --modification_probabilities: bed4 file containing modified CPG site probabilities [Required]
    -r, --regions: bed file containing HMM CDR and Transition Region Predictionsgit  [Required]
    -o, --output: Output file. Default is 'histogram.png' [Optional]
```

#### hmm_heatmaps.py
Generates two heatmaps to represent the emission and transition matrices from HMM CDR predictions. 
```
usage: hmm_heatmaps.py [-h]  
    -e, --emissionMatrix: Emmission matrix in .csv format from HMMCDRDetection.py [Required]
    -t, --transitionMatrix: Transition matrix in .csv format from HMMCDRDetection.py [Required]
    -o, --outputPrefix: Output prefix. Default is 'hmm_heatmap' [Optional]
```

#### create_track_pngs.py
Creates a png that has the methylation, strict CDRs, HMM CDRs, and H1L portion of the CenSat track. 
```
usage: create_track_pngs.py [-h]  
    -p, --mod_prob: Path to the 5mC probabilities bed file [Required]
    -s, --strict: Path to the Strict CDR bed file [Required]
    -v, --viterbi: Path to the Viterbi HMM bed file [Required]
    -r, --regions_file: Path to the satellite array regions bed file [Required]
    -o, --output_file: Output file path [Required]
```

