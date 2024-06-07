#!/bin/bash

#set -eux -o pipefail

show_help() {
    echo "usage: strict_score.sh [-h]
    -i, --input: Input bed4 containing fraction_modified information from modkit pileup. [required]
    -r, --hg002_merged_H1L: Bed Track of regions you want to find CDRs within. CDRs are identified by regions of aggregate depleted 5mC methylation. [required]
    -o, --output_prefix: A prefix for the output files. The output file will be named as <output_prefix>.strictScores.bed. [Required]
	--window_size: Set the size of the scoring windows. [default 510]
    --low_percent: The starting percentile for scoring windows. [default 1]
    --high_percent: Optional. The ending percentile for scoring windows. [default 20]
    --high-confidence: Cutoff for high confidence CDR score. Scores equal or above are high confidence. [default 15]
    --med-confidence: Cutoff for medium confidence CDR score. Scores above are medium confidence. [default 10]
    --low-confidence: Cutoff for low confidence CDR score. Scores above are low confidence. [default 3]"
}

# Default values
window_size=510
low_percent=1
high_percent=20
high_confidence=17
med_confidence=10
low_confidence=5

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h|--help) show_help; exit 0 ;;
        -i|--input) input="$2"; shift ;;
        -r|--hg002_merged_H1L) hg002_merged_H1L="$2"; shift ;;
        -o|--output_prefix) output_prefix="$2"; shift ;;
		--window_size) window_size="$2"; shift ;;
        --low_percent) low_percent="$2"; shift ;;
        --high_percent) high_percent="$2"; shift ;;
        --high-confidence) high_confidence="$2"; shift ;;
        --med-confidence) med_confidence="$2"; shift ;;
        --low-confidence) low_confidence="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; show_help; exit 1 ;;
    esac
    shift
done

# Check for required arguments
if [[ -z "$input" || -z "$hg002_merged_H1L" || -z "$output_prefix" ]]; then
    echo "Missing required arguments"
    show_help
    exit 1
fi

# Example processing of arguments
echo "Input file: $input"
echo "HG002 merged H1L: $hg002_merged_H1L"
echo "Output prefix: $output_prefix"
echo "Window size: $window_size"
echo "Low percentage: $low_percent"
echo "High percentage: $high_percent"
echo "High confidence score: $high_confidence"
echo "Medium confidence score: $med_confidence"
echo "Low confidence score: $low_confidence"

# make output folders
mkdir -p "windows"

# declare output variables
window_bed="windows/${output_prefix}.windows${window_size}.bed"
window_mean="windows/${output_prefix}.windows${window_size}.mean.bed"
cdr_scores="${output_prefix}.strictScores.bed"
strict_cdrs="${output_prefix}.strict.bed"

# generate the 1kb windows bed from the H1L censat annotations
awk -v win_size="$window_size" 'BEGIN { min = "unset"; max = 0 }
{
    if (min == "unset" || $2 < min) min = $2;
    if ($3 > max) max = $3;
}
END {
    # Assuming the chromosome is the same for all rows and is in the first column
    chrom = $1;
    for (i = min; i <= max; i += win_size) {
        window_start = i;
        window_end = (i + win_size - 1 > max) ? max : i + win_size - 1;
        print chrom "\t" window_start "\t" window_end;
    }
}' $input > $window_bed

# generate bed
bedtools intersect -a $window_bed -b $hg002_merged_H1L | \
	sort -k 1,1 -k2,2n - | \
	bedtools map -a - -b $input -c 4 -o mean | \
	awk -F'\t' '$4 != "." {print}' - | \
	sort -k 1,1 -k2,2n - > $window_mean

# initialize cdr scoring bedgraph
awk 'OFS="\t" {print $1, $2, $3, 0}' $window_bed > $cdr_scores

# Store scoring bed file in an associative array
declare -A windows
while read -r chrom start end score; do
	windows["$chrom,$start,$end"]=$score
done < $cdr_scores

for ((percent=high_percent; percent>=low_percent; percent--)); do
	# generate the thresholds for the strict CDR/transitions from the percentages
	threshold=$(awk '{print $4}' $window_mean | sort -n | \
		awk -v perc=$percent 'BEGIN{line=-1} {all[NR]=$1} END{line=int((perc/100.0)*NR); if(line<1)line=1; print all[line]}')

	# Loop through each line of the methylation bed file
	while IFS=$'\t' read -r chrom start end methylation; do
		# Check if the average methylation is below the threshold
		if (( $(echo "$methylation < $threshold" | bc -l) )); then
			windows["$chrom,$start,$end"]=$((windows["$chrom,$start,$end"] + 1))
		fi
	done < $window_mean
done

# Write updated scoring bed file
for key in "${!windows[@]}"; do
    # Replace commas with tabs in the key
    formatted_key="${key//,/$'\t'}"
    # Echo the formatted key-value pair
    echo -e "$formatted_key\t${windows[$key]}"
done > tmpfile && mv tmpfile ${cdr_scores}

# sort output :)
sort -k 1,1 -k2,2n -o ${cdr_scores} ${cdr_scores}

# Find the Strict High Confidence CDRs (merge nearby high confidence windows with this)
awk -v high_conf=$high_confidence 'BEGIN {OFS=FS="\t"} $4 >= high_conf {print $1, $2, $3, $4}' ${cdr_scores} | \
	sort -k 1,1 -k2,2n - | \
	bedtools merge -d $(echo "4 * $window_size + 1" | bc) -c 4 -o 'mean' -i - | \ 
	awk 'BEGIN {OFS=FS="\t"} $3 - $2 >= 1750 {print $1, $2, $3, "strict_CDR", 0, ".", $2, $3, "0,0,255"}' | \
	sort -k 1,1 -k2,2n - | \
	awk '!seen[$0]++' > temp_cdrs.bed

# Find the Strict Medium Confidence CDRs (only merge adjacent windows)
awk -v med_conf=$med_confidence 'BEGIN {OFS=FS="\t"} $4 > med_conf {print $1, $2, $3, $4}' ${cdr_scores} | \
	sort -k 1,1 -k2,2n - | \
	bedtools merge -d 1 -c 4 -o 'mean' -i - | \
	bedtools subtract -a - -b temp_cdrs.bed | \
	awk 'BEGIN {OFS=FS="\t"} $3 - $2 >= 1750 {print $1, $2, $3, "strict_CDR", 0, ".", $2, $3, "0,0,255"}' | \
	sort -k 1,1 -k2,2n - | \
	awk '!seen[$0]++' >> temp_cdrs.bed
	
# Find the Strict Transitions
awk -v low_conf=$low_confidence 'BEGIN {OFS=FS="\t"} $4 > low_conf {print $1, $2, $3, $4}' ${cdr_scores} | \
	sort -k 1,1 -k2,2n - | \
	bedtools merge -d $(echo "2 * $window_size + 1" | bc) -c 4 -o 'mean' -i - | \
	bedtools intersect -a - -b temp_cdrs.bed -wa | \
	bedtools subtract -a - -b temp_cdrs.bed | \
	awk 'BEGIN {OFS=FS="\t"} {print $1, $2, $3, "strict_Transition", 0, ".", $2, $3, "173,216,230"}' | \
	sort -k 1,1 -k2,2n - | \
	awk '!seen[$0]++' > temp_transitions.bed

# generate output
cat temp_cdrs.bed temp_transitions.bed | \
	sort -k1,1 -k2,2n -o ${strict_cdrs} 

rm temp_cdrs.bed &&
    rm temp_transitions.bed

echo "Strict CDR File Written to: ${strict_cdrs}"