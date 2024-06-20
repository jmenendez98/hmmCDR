#!/bin/bash

set -eux -o pipefail

# Initialize variables
file=""
hg002_merged_H1L=""
prefix=""
percent=10
transition_percent=20
window_size=1020
min_size=3
merge_distance=3

OPTIONS=i:,r:,o:,p:,t:,w:,m:,d:
LONGOPTS=input:,hg002_merged_H1L:,output_prefix:,percent:,transition_percent:,window_size:,min_size:,merge_distance:

# Check if the correct number of arguments are provided
if [ "$#" -lt 1 ]; then
    echo "Error: Missing required arguments. Use -h for help."
    exit 1
fi

# Parse command-line options
PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ $? -ne 0 ]]; then
    exit 2
fi

# Using eval to handle options correctly
eval set -- "$PARSED"

# Extract options and their arguments into variables
while true; do
    case "$1" in
        -i|--input)
            file="$2"
            shift 2
            ;;
        -r|--hg002_merged_H1L)
            hg002_merged_H1L="$2"
            shift 2
            ;;
        -o|--output_prefix)
            prefix=$(basename "$2")
            shift 2
            ;;
        -p|--percent)
            percent="$2"
            shift 2
            ;;
        -t|--transition_percent)
            transition_percent="$2"
            shift 2
            ;;
        -w|--window_size)
            window_size="$2"
            shift 2
            ;;
        -m|--min_size)
            min_size="$2"
            shift 2
            ;;
		-d|--merge_distance)
            merge_distance="$2"
            shift 2
            ;;	
        -h|--help)
            echo "Usage: strictCDRDetection.sh -i <input_file> -r <hg002_merged_H1L_file> -o <output_prefix> -p <percentage> -t <transition_percentage> -w <window_size> -m <min_size>"
            echo "Options:"
            echo "  -i, --input: Required. Modkit pileup bedgraph of fraction modified."
            echo "  -r, --hg002_merged_H1L: Required. The file containing the H1L(active-Alpha Satellite) regions."
            echo "  -o, --output_prefix: Required. A prefix for the output files. The output files will be named as <output_prefix>.strictCDR.bed and <output_prefix>.strictTransitions.bed."
            echo "  -p, --percent: Optional. The percentage threshold for the CDRs. Default is 10."
            echo "  -t, --transition_percent: Optional. The transition percentage threshold. Default is 20."
            echo "  -w, --window_size: Optional. The window size. Default is 1020."
            echo "  -m, --min_size: Optional. The minimum amound of windows flagged as a strict CDRs to keep a CDR. Minimum CDR Size is this * window_size. Default is 3."
			echo "  -d, --merge_distance: Optional. The distance to merge nearby CDR entries(before filtering by size). Distance is given in number of windows. Default is 2."
            exit 0
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Unexpected option: $1"
            exit 3
            ;;
    esac
done

# Check if required options are provided
if [[ -z "$file" || -z "$hg002_merged_H1L" || -z "$prefix" ]]; then
	echo "Missing required arguments"
	exit 1
fi

# make output folders
mkdir -p "windows"

# declare variables
window_bed="windows/${prefix}.windows1000.bed"
window_mean="windows/${prefix}.windows1000.mean.bed"
strict_cdrs="${prefix}.strictCDRs.bed"

# 1
awk -v window_size="$window_size" '
BEGIN { min = "unset"; max = 0 }
{
	if (min == "unset" || $2 < min) min = $2;
	if ($3 > max) max = $3;
}
END {
	# Assuming the chromosome is the same for all rows and is in the first column
	chrom = $1;
	for (i = min; i <= max; i += window_size) {
		window_start = i;
		window_end = (i + window_size - 1 > max) ? max : i + window_size - 1;
		print chrom "\t" window_start "\t" window_end;
	}
}' $file > $window_bed

# 2 + 3
bedtools intersect -a $window_bed -b $hg002_merged_H1L | \
	sort -k 1,1 -k2,2n - | \
	bedtools map -a - -b $file -c 4 -o mean | \
	awk -F'\t' '$4 != "." {print}' - | \
	sort -k 1,1 -k2,2n - > $window_mean

# Reset current thresholds for each file! 
current_percent=$percent
current_transition_percent=$transition_percent

# 4
cdr_threshold=$(awk '{print $4}' $window_mean | sort -n | \
	awk -v perc=$current_percent 'BEGIN{line=-1} {all[NR]=$1} END{line=int((perc/100.0)*NR); if(line<1)line=1; print all[line]}')

cdr_transition_threshold=$(awk '{print $4}' $window_mean | sort -n | \
	awk -v perc=$current_transition_percent 'BEGIN{line=-1} {all[NR]=$1} END{line=int((perc/100.0)*NR); if(line<1)line=1; print all[line]}')

merge_dist=$(echo "$merge_distance * $window_size + 1" | bc)
min_cdr_size=$(echo "$min_size * $window_size + 1" | bc)

# 5
awk -v cdr_thresh=$cdr_threshold '$4 <= cdr_thresh' $window_mean | \
	bedtools merge -d 3 -i - | \
	bedtools intersect -a - -b $hg002_merged_H1L -f 1.0 | \
	awk -v min=$min_cdr_size 'BEGIN {FS=OFS="\t"} {if ($3-$2 > min) {print $1,$2,$3}}' - | \
	bedtools merge -d $merge_dist -i - | \
    awk 'BEGIN {OFS=FS="\t"} $3 - $2 >= 1750 {print $1, $2, $3, "strict_CDR", 0, ".", $2, $3, "0,0,255"}' | \
	awk '!seen[$0]++' > temp_cdrs.bed

awk -v trans_thresh=$cdr_transition_threshold 'BEGIN {FS=OFS="\t"} $4 <= trans_thresh' $window_mean | \
	bedtools merge -d 3 -i - | \
	bedtools intersect -a - -b $hg002_merged_H1L -f 1.0 | \
	bedtools intersect -a - -b temp_cdrs.bed -wa | \
	bedtools subtract -a - -b temp_cdrs.bed | \
    awk 'BEGIN {OFS=FS="\t"} {print $1, $2, $3, "strict_Transition", 0, ".", $2, $3, "173,216,230"}' | \
	awk '!seen[$0]++' > temp_transitions.bed

cat temp_cdrs.bed temp_transitions.bed | \
    sort -k 1,1 -k2,2n -o ${strict_cdrs}

rm temp_cdrs.bed temp_transitions.bed

echo "Wrote to: ${strict_cdrs}"