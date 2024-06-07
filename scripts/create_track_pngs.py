import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


# Function to parse bed file
def parse_pileup_file(bed_file):
    data = []
    with open(bed_file, 'r') as f:
        for line in f:
            chrom, start, end, prob = line.strip().split('\t')
            data.append((chrom, int(start), int(end), float(prob)))
    return data

# Function to parse pileup bed file
def parse_strict_file(bed_file):
    data = []
    with open(bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 4:
                chrom, start, end, strict_type = fields[:4]
                data.append((chrom, int(start), int(end), strict_type))
    return data

# Function to parse bed file
def parse_viterbi_file(bed_file):
    data = []
    with open(bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 4:
                chrom, start, end, cdr_type = fields[:4]
                data.append((chrom, int(start), int(end), cdr_type))
    return data

# Function to parse regions file
def parse_regions_file(regions_file):
    data = []
    with open(regions_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 4:
                chrom, start, end, array_type = fields[:4]
                data.append((chrom, int(start), int(end), array_type))
    return data


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Generate IGV-like visualization of 5mC probabilities and satellite array regions.')
    parser.add_argument('-p', '--mod_prob', help='Path to the 5mC probabilities bed file')
    parser.add_argument('-s', '--strict', help='Path to the Strict CDR bed file')
    parser.add_argument('-v', '--viterbi', help='Path to the Viterbi HMM bed file')
    parser.add_argument('-r', '--regions_file', help='Path to the satellite array regions bed file')
    parser.add_argument('-o', '--output_file', help='Output file path')
    args = parser.parse_args()

    # Parse Each File
    pileup_data = parse_pileup_file(args.mod_prob)
    strict_data = parse_strict_file(args.strict)
    viterbi_data = parse_viterbi_file(args.viterbi)
    regions_data = parse_regions_file(args.regions_file)

    # Plot settings
    fig, ax = plt.subplots(figsize=(10, 2.25))

    #########################
    ### ADD ITEMS TO PLOT ###
    #########################

    # Plot the locations of the H1Ls from the CenSat annotation track
    xmin = np.inf # min and max determined by CenSat annotations
    xmax = -np.inf
    for chrom, start, end, array_type in regions_data:
        if chrom == chromosome:
            if 'H1L' in array_type:
                ax.plot([start, end], 
                        [(panel_step*6), (panel_step*6)], 
                        color='red', 
                        linewidth=2)
                if start < xmin:
                    xmin = start
                if end > xmax:
                    xmax = end

    print(chromosome)
    print(xmin, xmax)

    panel_step = 0.25

    # Plot the bedgraph of CpG probabilities
    for chrom, start, end, prob in pileup_data:
        width = end - start  # Calculate the width of the bar
        rect = Rectangle(xy=(start, panel_step*3), 
                         width=width, 
                         height=((prob*0.01)*(panel_step*2)), 
                         color='blue',
                         alpha=0.025)
        ax.add_patch(rect)
        chromosome = chrom

    # Plot the strictCDRs from the bed file of strict CDR Predictions
    strict_color = {"strict_CDR": "orange",
                    "strict_Transition": "navajowhite"}
    for chrom, start, end, strict_type in strict_data:
        if chrom == chromosome:
            ax.plot([start, end], 
                    [panel_step, panel_step], 
                    color=strict_color[strict_type], 
                    linewidth=2)          


    # Plot the CDRs from the Viterbi HMM CDR Predictions
    cdr_color = {"CDR": "maroon",
                 "small_CDR": "darksalmon",
                 "CDR_Transition": "lightsalmon"}
    for chrom, start, end, cdr_type in viterbi_data:
        if chrom == chromosome:
            ax.plot([start, end], 
                    [(panel_step*2), (panel_step*2)], 
                    color=cdr_color[cdr_type], 
                    linewidth=2)  



    #########################
    ### SET PLOT SETTINGS ###
    #########################

    # set x-axis limits to 1.1 times the length of the array
    #middle = (xmin+xmax)/2
    #half_xwidth = (abs(xmax - xmin)/2)*1.05
    #ax.set_xlim(middle-half_xwidth, middle+half_xwidth)
    ax.set_xlim(xmin, xmax)
    # print("H1L Limits: ", xmin, xmax)

    # Set x-axis tick locations and labels
    interval_size = (xmax - xmin) // 10
    tick_locations = np.arange(xmin, xmax, step=interval_size)
    tick_labels = [str(pos) for pos in tick_locations]  # Convert positions to string labels
    plt.xticks(tick_locations, tick_labels, rotation=45)
        
    ax.set_ylim(0, panel_step*(panel_step*7))

    # Set labels and title
    #ax.set_xlabel('Genomic Position')
    #ax.set_ylabel('5mC Probability')
    ax.set_title(chromosome)

    #ax.yaxis.set_ticks([])  # Remove Y axis ticks
    #ax.yaxis.set_ticklabels([])  # Remove Y axis tick labels
    ax.yaxis.set_ticks([panel_step, panel_step*2, panel_step*4,panel_step*6])  # Place labels of Y 
    ax.yaxis.set_ticklabels(["Strict CDR", "Viterbi HMM CDR", "Modkit Fraction Modified", "CenSat Live Array"])  # Label Tracks 

    # Remove top, left, and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # remove y axis ticks
    ax.tick_params(axis='y', left=False)

    # Show plot
    plt.tight_layout()
    plt.savefig(args.output_file, dpi=2400)  # Save plot as PNG image
    plt.show()

if __name__ == "__main__":
    main()
